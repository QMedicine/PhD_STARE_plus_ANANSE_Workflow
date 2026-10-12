#!/bin/bash
#SBATCH --job-name=find_TADs
#SBATCH --array=0-17
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicFindTADs_logs_reps_merged/find_TADs_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicFindTADs_logs_reps_merged/find_TADs_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

module load miniconda3
conda activate hicexplorer

# Directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files/hicFindTADs_TADs"
mkdir -p ${OUTPUT_DIR}

# List all corrected .h5 files
FILES=($(ls ${INPUT_DIR}/*.h5 | sort))

# Select the file for this task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${FILE}" .h5)

# Extract the resolution from the file name
RESOLUTION=$(echo ${BASENAME} | awk -F'_' '{print $3}')
echo "Processing file: ${FILE}"
echo "Extracted resolution: ${RESOLUTION}"

# Calculate minDepth and maxDepth based on resolution
MIN_DEPTH=$((RESOLUTION * 3))  # Minimum depth = 3x resolution
MAX_DEPTH=$((RESOLUTION * 10))  # Maximum depth = 10x resolution
echo "Calculated minDepth: ${MIN_DEPTH}"
echo "Calculated maxDepth: ${MAX_DEPTH}"

# Run hicFindTADs
hicFindTADs \
    --matrix "${FILE}" \
    --outPrefix "${OUTPUT_DIR}/${BASENAME}" \
    --minDepth ${MIN_DEPTH} \
    --maxDepth ${MAX_DEPTH} \
    --step ${RESOLUTION} \
    --thresholdComparisons 0.05 \
    --delta 0.01 \
    --correctForMultipleTesting fdr \
    -p 4  # Use 4 processors for parallel computation

conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicFindTADs_replicates_merged.sh