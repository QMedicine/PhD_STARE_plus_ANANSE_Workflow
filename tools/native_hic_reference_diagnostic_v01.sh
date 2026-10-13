#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_reference_check
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" python - <<'PY'
from pathlib import Path
import pysam,subprocess,json
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
bam=root/'runs/hic_single_library_pilot_v36/results/raw/HiC/O_LTHSC_Rep1/O_LTHSC_Rep1.name_sorted.bam'
with pysam.AlignmentFile(bam,'rb') as f:refs=list(zip(f.references,f.lengths))
fasta=Path('/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa')
names=[]
with fasta.open() as f:
 for line in f:
  if line.startswith('>'):names.append(line.split()[0][1:])
result=dict(bam_reference_count=len(refs),bam_references=refs,fasta_reference_count=len(names),fasta_references=names,bowtie_summary=subprocess.check_output(['bowtie2-inspect','-s','/vast/projects/Immune_Aging/Reference_genomes/mm39/Bowtie2/mm39'],text=True))
out=root/'evidence/hic_reference_diagnostic_v01';out.mkdir(exist_ok=True);(out/'reference_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
PY
