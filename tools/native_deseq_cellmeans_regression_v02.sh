#!/usr/bin/env bash
#SBATCH --job-name=phd_deseq_cellmeans
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
run="$TASK_ROOT/runs/deseq_cellmeans_v02_$SLURM_JOB_ID";mkdir -p "$run";cd "$run"
image=/vast/projects/Immune_Aging/apptainer_singularity_cache/depot.galaxyproject.org-singularity-bioconductor-deseq2-1.34.0--r41hc247a5b_3.img
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" Rscript "$TASK_ROOT/source_v07/bin/deseq_historical.R" \
 "$TASK_ROOT/runs/rna_recount_v06/results/rna/tpm/counts.counts_filtered.tsv" \
 /vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/samplesheet_Old_and_Young_HSCs_v01.tsv \
 Old_HSC Young_HSC forward.deseq2.tsv reverse.deseq2.tsv 1 'Old_HSC,Young_HSC'
image=/vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img
expected=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" python "$TASK_ROOT/tools/compare_numeric_tables.py" \
 --expected "$expected/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv" --observed forward.deseq2.tsv --output forward_comparison.json
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" python "$TASK_ROOT/tools/compare_numeric_tables.py" \
 --expected "$expected/HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv" --observed reverse.deseq2.tsv --output reverse_comparison.json
cat forward_comparison.json reverse_comparison.json
