#!/usr/bin/env bash
#SBATCH --job-name=phd_rna_raw_regression_v02
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/rna_fastq_regression_v02";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$TASK_ROOT/tools/native_rna_fastq_regression_v02.py" \
 --run "$TASK_ROOT/runs/raw_RNA_integration_v14" --reference-run "$TASK_ROOT/runs/rna_recount_v08" \
 --reference-counts /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt \
 --trimmed-dir /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed --output rna_fastq_regression.json
