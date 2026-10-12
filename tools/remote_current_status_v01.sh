set -euo pipefail
squeue -u stacpoole.q -o '%.18i %.35j %.8T %.10M %R'
for file in /vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/logs/rna_recount_*.err /vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/logs/rna_recount_*.out; do
 printf '\nLOG %s\n' "$file";tail -n 15 "$file"
done
script_directory="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Hannah's_Hi-C_fastq_alignment_QC_processing_scripts"
ls -la "$script_directory"
