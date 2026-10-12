#!/usr/bin/env bash
#SBATCH --job-name=phd_rna_regression
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/rna_numerical_regression_v01";mkdir -p "$run";cd "$run"
image=/vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img
observed="$TASK_ROOT/runs/rna_recount_v06/results/rna"
expected=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results
compare(){
 apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" \
  python "$TASK_ROOT/tools/compare_numeric_tables.py" --expected "$1" --observed "$2" --output "$3"
}
compare /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv "$observed/recount/selected_counts.csv" cleaned_counts.json
compare "$expected/HSC_50kbp_window_H13CORE_motifs/01_counts/02_tpm/counts.tpm.tsv" "$observed/tpm/counts.tpm.tsv" unfiltered_tpm.json
compare "$expected/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv" "$observed/de/Old_HSC_to_Young_HSC/forward.deseq2.tsv" forward_de.json
compare "$expected/HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv" "$observed/de/Old_HSC_to_Young_HSC/reverse.deseq2.tsv" reverse_de.json
cat cleaned_counts.json unfiltered_tpm.json forward_de.json reverse_de.json
