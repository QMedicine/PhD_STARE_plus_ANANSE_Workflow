#!/bin/bash
#SBATCH --job-name=correct_h5
#SBATCH --array=0-17
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicCorrectMatrix_logs_reps_merged/hicCorrectMatrix_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicCorrectMatrix_logs_reps_merged/hicCorrectMatrix_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# Input and output directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files"
THRESHOLDS_TSV="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_diag_plots/thresholds_table_replicates_merged.tsv"
mkdir -p ${OUTPUT_DIR}

# Get the list of files
FILES=($(ls ${INPUT_DIR}/*.h5))

# Select the current file based on the Slurm array task ID
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${FILE}" .h5)

# Extract thresholds from the TSV file
    THRESHOLDS=$(awk -v filename="${BASENAME}.h5" '$1 == filename {print $3, $4}' ${THRESHOLDS_TSV})
    LOWER_THRESH=$(echo ${THRESHOLDS} | cut -d' ' -f1)
    UPPER_THRESH=$(echo ${THRESHOLDS} | cut -d' ' -f2)


# Echo checks for verification
echo "Processing file: ${FILE}"
echo "Lower threshold: ${LOWER_THRESH}"
echo "Upper threshold: ${UPPER_THRESH}"

# Correct the matrix
hicCorrectMatrix correct \
    --matrix "${FILE}" \
    --correctionMethod ICE \
    --filterThreshold ${LOWER_THRESH} ${UPPER_THRESH} \
    --outFileName "${OUTPUT_DIR}/${BASENAME}_hicCorrectMatrix.h5"

# Deactivate conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicCorrectMatrix_correct_replicates_merged.sh