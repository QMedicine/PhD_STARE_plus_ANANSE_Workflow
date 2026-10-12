#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_qualify_v03
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
image="$TASK_ROOT/containers/hicexplorer_v01.sif"
run="$TASK_ROOT/runs/hic_qualification_v03_$SLURM_JOB_ID"
mkdir -p "$run";cd "$run"
sha256sum "$image" > image.sha256
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" python - <<'PY'
import cooler,pandas as pd,json
from pathlib import Path
chroms=[str(x) for x in range(1,20)]+['X']
bins=pd.DataFrame([dict(chrom=c,start=i*100000,end=(i+1)*100000) for c in chroms for i in range(2)])
pixels=pd.DataFrame([dict(bin1_id=2*i,bin2_id=2*i+1,count=5) for i in range(20)])
cooler.create_cooler('tiny.cool',bins,pixels)
source='/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Coolers/OandY_LTHSC_MMP_Merged_Cool_files/O_LTHSC_Rep1_100000_MMP_combined.cool'
c=cooler.Cooler(source)
Path('historical_cooler_metadata.json').write_text(json.dumps(dict(info=c.info,chromsizes={k:int(v) for k,v in c.chromsizes.items()}),indent=2,default=lambda v:v.item())+'\n')
PY
apptainer exec --cleanenv --no-home "$image" hicConvertFormat --matrices tiny.cool --inputFormat cool --outputFormat h5 --outFileName tiny.h5
apptainer exec --cleanenv --no-home "$image" hicConvertFormat --matrices tiny.h5 --inputFormat h5 --outputFormat ginteractions --outFileName contacts_all
ls -l
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$image" python "$TASK_ROOT/source_v03/bin/split_contacts.py" --input contacts_all.tsv --outdir contacts
apptainer exec --cleanenv --no-home "$image" python - <<'PY'
import gzip
from pathlib import Path
files=list(Path('contacts').glob('*.gz'))
assert len(files)==20
for p in files:
 with gzip.open(p,'rt') as f:
  rows=[line.strip().split('\t') for line in f];assert len(rows)==1;assert rows[0][:2]==['0','100000'];assert float(rows[0][2])==5.0
print('FUNCTIONAL_CONTACT_TEST_PASS')
PY
printf 'HIC_QUALIFICATION_COMPLETE\n'
