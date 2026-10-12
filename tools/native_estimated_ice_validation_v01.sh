#!/usr/bin/env bash
#SBATCH --job-name=phd_estimated_ice_validate
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/estimated_ice_validation_v01";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/hicexplorer_v01.sif" python - <<'PY'
from pathlib import Path
import hashlib,json,numpy as np
from hicmatrix import HiCMatrix
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/runs/cooler_estimated_ice_v15/results/estimated_ice')
records=[]
for path in sorted(root.glob('*/*.estimated_ICE.h5')):
 obj=HiCMatrix.hiCMatrix(str(path));data=obj.matrix.data;delta=(obj.matrix-obj.matrix.T).data
 if not len(data) or not np.isfinite(data).all() or np.min(data)<0 or (len(delta) and np.max(np.abs(delta))>1e-8):raise ValueError('Estimated matrix is empty, nonfinite, negative or asymmetric')
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 provenance=json.loads((path.parent/'hic_estimated_provenance.json').read_text())
 if provenance['status']!='ESTIMATED_PRE_BLACKLIST_NOT_HISTORICALLY_AUTHENTICATED':raise ValueError('Missing estimated provenance label')
 records.append(dict(cell_state=provenance['cell_state'],n_bins=len(obj.cut_intervals),nonzero_entries=int(obj.matrix.nnz),min_value=float(np.min(data)),max_value=float(np.max(data)),all_finite=True,symmetric_within_1e_8=True,sha256=h.hexdigest(),bytes=path.stat().st_size,status=provenance['status']))
if len(records)!=2:raise ValueError('Expected two study estimated matrices')
Path('estimated_ice_validation.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
PY
