#!/usr/bin/env bash
#SBATCH --job-name=phd_rna_fastq_regression
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/rna_fastq_regression_v01";mkdir -p "$run";cd "$run"
image="$TASK_ROOT/containers/raw_tools_v06.sif"
source="$TASK_ROOT/runs/raw_RNA_pilot_v13/results/raw/RNA/SRR13192285_Young_HSC_Rep1"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" featureCounts -T 8 \
 -a /vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf --extraAttributes gene_name,gene_biotype \
 -o recounted.featureCounts.txt "$source/SRR13192285_Young_HSC_Rep1_Aligned.sortedByCoord.out.bam"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" python - <<'PY'
from pathlib import Path
import gzip,hashlib,json
import pandas as pd
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
pilot=root/'runs/raw_RNA_pilot_v13/results/raw/RNA/SRR13192285_Young_HSC_Rep1'
a=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed/SRR13192285_Young_HSC_Rep1_trimmed.fastq.gz');b=pilot/'SRR13192285_Young_HSC_Rep1.trimmed_R1.fastq.gz'
def decoded_sha(p):
 h=hashlib.sha256();n=0
 with gzip.open(p,'rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block);n+=len(block)
 return h.hexdigest(),n
sha_a,size_a=decoded_sha(a);sha_b,size_b=decoded_sha(b)
expected=pd.read_csv('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt',sep='\t',comment='#',index_col=0)
expected.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in expected.columns]
observed=pd.read_csv('recounted.featureCounts.txt',sep='\t',comment='#',index_col=0)
if set(expected.index)!=set(observed.index):raise ValueError('featureCounts gene identities differ')
x=expected['SRR13192285_Young_HSC_Rep1'];y=observed.iloc[:,-1].loc[expected.index]
record=dict(sample='SRR13192285_Young_HSC_Rep1',trimmed_sequence_bytes_identical=sha_a==sha_b,expected_decoded_sha256=sha_a,observed_decoded_sha256=sha_b,expected_decoded_bytes=size_a,observed_decoded_bytes=size_b,n_featurecounts_genes=len(x),different_gene_count_values=int((x!=y).sum()),counts_exact=bool((x==y).all()))
Path('rna_fastq_regression.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if not record['trimmed_sequence_bytes_identical'] or not record['counts_exact']:raise SystemExit('RNA raw regression differs')
PY
