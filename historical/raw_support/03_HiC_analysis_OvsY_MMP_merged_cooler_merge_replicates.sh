#!/bin/bash
#SBATCH --job-name=merge_cool_replicates
#SBATCH --array=0-17
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/merged_replicates_logs/merge_cool_replicates_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/merged_replicates_logs/merge_cool_replicates_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Activate conda environment
module load miniconda3
conda activate cooler-env

# Directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_1000bp_Cool_files"
OUTPUT_DIR="${INPUT_DIR}/merged_replicates_cool_files"
mkdir -p "${OUTPUT_DIR}"

# Extract unique condition-resolution combinations
CONDITIONS_RESOLUTIONS=($(ls ${INPUT_DIR}/*_MMP_combined_BL_removed.cool | xargs -n 1 basename | awk -F'_' '{print $1"_"$4}' | sort -u))

# Select the current condition-resolution based on the array index
CONDITION_RES=${CONDITIONS_RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}
CONDITION=$(echo ${CONDITION_RES} | awk -F'_' '{print $1}')
RESOLUTION=$(echo ${CONDITION_RES} | awk -F'_' '{print $2}')

echo "Merging replicates for condition: ${CONDITION}, resolution: ${RESOLUTION}"

# Find files to merge
FILES_TO_MERGE=($(ls ${INPUT_DIR}/${CONDITION}_LTHSC_Rep*_${RESOLUTION}_MMP_combined_BL_removed.cool 2>/dev/null))

if [[ ${#FILES_TO_MERGE[@]} -eq 3 ]]; then
    echo "  Files to merge:"
    for FILE in "${FILES_TO_MERGE[@]}"; do
        echo "    ${FILE}"
    done

    # Define output file
    OUTPUT_FILE="${OUTPUT_DIR}/${CONDITION}_LTHSC_${RESOLUTION}_MMP_combined_replicates_merged_BL_removed.cool"

    # Merge files
    cooler merge "${OUTPUT_FILE}" "${FILES_TO_MERGE[@]}"

    echo "Merged file saved to: ${OUTPUT_FILE}"
else
    echo "  Error: Expected 3 replicates, found ${#FILES_TO_MERGE[@]}. Skipping."
fi

conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_cooler_merge_replicates.sh