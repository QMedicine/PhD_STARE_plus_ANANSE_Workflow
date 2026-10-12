#!/usr/bin/env bash
#SBATCH --job-name=phd_atac_header_diagnostic
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
run="$TASK_ROOT/evidence/atac_header_diagnostic_v01";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/raw_tools_v06.sif" python - <<'PY'
import hashlib,json,subprocess
from pathlib import Path
name='SRR13186609_Young_HSC_B1_Rep1';root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
paths={'original':Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered')/(name+'_filtered.bam'),'rebuilt':root/'runs/raw_ATAC_pilot_v15/results/raw/ATAC'/name/(name+'.filtered.bam')}
record={}
for label,path in paths.items():
 header=subprocess.check_output(['samtools','view','-H',str(path)],text=True)
 Path(label+'.header.txt').write_text(header)
 reader=subprocess.Popen(['samtools','view',str(path)],stdout=subprocess.PIPE,text=True)
 rows=[]
 for i,line in enumerate(reader.stdout):
  values=line.split('\t');rows.append(dict(qname=values[0],flag=values[1],rname=values[2],pos=values[3],mapq=values[4],cigar=values[5],sequence_sha256=hashlib.sha256(values[9].encode()).hexdigest(),quality_sha256=hashlib.sha256(values[10].encode()).hexdigest()))
  if i==9:break
 reader.stdout.close();reader.terminate();reader.wait()
 record[label]=dict(program_headers=[line for line in header.splitlines() if line.startswith('@PG')],first10_records=rows)
Path('atac_header_diagnostic.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({label:value['program_headers'] for label,value in record.items()},indent=2))
PY
