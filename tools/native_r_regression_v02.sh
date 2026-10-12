#!/usr/bin/env bash
#SBATCH --job-name=phd_r_regression_v02
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
image="$TASK_ROOT/containers/r_plot_v01.sif"
sha256sum --check "$TASK_ROOT/evidence/r_plot_v01.sif.sha256"
run="$TASK_ROOT/runs/r_regression_v02_$SLURM_JOB_ID"
mkdir -p "$run";cd "$run"
printf 'JOB=%s HOST=%s ARCH=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$(uname -m)"
an=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM
rn=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results
for direction in Old2Young Young2Old; do
 if [[ "$direction" == Old2Young ]]; then
  source=Old_HSC;target=Young_HSC;stem=STARE_OldHSC_to_YoungHSC_TF_influence
  de="$rn/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv"
 else
  source=Young_HSC;target=Old_HSC;stem=STARE_YoungHSC_to_OldHSC_TF_influence
  de="$rn/HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv"
 fi
 mkdir -p "$direction"
 apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" Rscript "$TASK_ROOT/source_v02/bin/plot_top20.R" \
  "$TASK_ROOT/source_v02/historical/scripts/09_plots/GRN_plotting_${direction}.Rmd" \
  "$an/TF_Influence_$direction/$stem.txt" "$an/TF_Influence_$direction/${stem}_diffnetwork.tsv" "$de" "$source" "$target" "$direction/plots"
done
sha256sum Old2Young/plots/*.pdf Young2Old/plots/*.pdf > pdf_checksums.sha256
printf 'R_REGRESSION_COMPLETE\n'
