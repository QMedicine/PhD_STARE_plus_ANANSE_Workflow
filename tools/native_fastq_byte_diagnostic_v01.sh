#!/usr/bin/env bash
#SBATCH --job-name=phd_fastq_byte_diagnostic
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/fastq_byte_diagnostic_v01";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" python - <<'PY'
import hashlib,json,numpy as np
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/inputs')
paths=[root/'public_fastq/SRR13186609.fastq.gz.partial',root/'public_fastq_ftp_v02/SRR13186609.fastq.gz.partial']
offset=0;different=0;first=[];chunks=[];hashes=[hashlib.md5(),hashlib.md5()]
with paths[0].open('rb') as a,paths[1].open('rb') as b:
 while True:
  left=a.read(4*1024*1024);right=b.read(4*1024*1024)
  if not left and not right:break
  if len(left)!=len(right):raise ValueError('Failed transfers have different sizes')
  hashes[0].update(left);hashes[1].update(right)
  indices=np.flatnonzero(np.frombuffer(left,dtype=np.uint8)!=np.frombuffer(right,dtype=np.uint8))
  different+=len(indices)
  if len(indices):
   chunks.append(dict(start=offset,end=offset+len(left),different_bytes=int(len(indices))))
   first.extend(int(i)+offset for i in indices[:max(0,20-len(first))])
  offset+=len(left)
record=dict(paths=[str(p) for p in paths],bytes=offset,md5_in_independent_container=[h.hexdigest() for h in hashes],different_bytes=int(different),first_different_offsets=first,different_chunks=chunks)
Path('byte_diagnostic.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
PY
