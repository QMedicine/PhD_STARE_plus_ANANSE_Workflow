#!/usr/bin/env bash
#SBATCH --job-name=phd_selected_plot_v04
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
SOURCE_DIR="$TASK_ROOT/source_v04"
run="$TASK_ROOT/runs/selected_plot_v04_$SLURM_JOB_ID"
mkdir -p "$run"; cd "$run"
an=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM
rn=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results
failed=0
for direction in Old2Young Young2Old; do
 if [[ "$direction" == Old2Young ]]; then
  source=Old_HSC;target=Young_HSC;stem=STARE_OldHSC_to_YoungHSC_TF_influence
  de="$rn/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv"
 else
  source=Young_HSC;target=Old_HSC;stem=STARE_YoungHSC_to_OldHSC_TF_influence
  de="$rn/HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv"
 fi
 mkdir -p "$direction"
 apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/r_plot_v01.sif" \
   Rscript "$SOURCE_DIR/bin/plot_top20.R" "$SOURCE_DIR/bin/plot_branching_template.R" \
   "$an/TF_Influence_$direction/$stem.txt" "$an/TF_Influence_$direction/${stem}_diffnetwork.tsv" "$de" "$source" "$target" "$direction/plots" || failed=1
 mv versions.txt "$direction/versions.txt"
done
printf 'SELECTED_PLOT_TEST_EXIT=%s\n' "$failed"
exit "$failed"
