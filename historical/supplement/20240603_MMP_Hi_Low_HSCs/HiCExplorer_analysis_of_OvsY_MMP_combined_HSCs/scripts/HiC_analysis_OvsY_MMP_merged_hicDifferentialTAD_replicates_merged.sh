#!/bin/bash
#SBATCH --job-name=hicDifferentialTAD
#SBATCH --array=0-8
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicDifferentialTAD_logs/hicDifferentialTAD_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicDifferentialTAD_logs/hicDifferentialTAD_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# Activate the conda environment for HiCExplorer
module load miniconda3
conda activate hicexplorer

# Input directories
MATRIX_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files"
TAD_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files/hicFindTADs_TADs"

# Output directory
OUTPUT_DIR="${TAD_DIR}/hicDifferentialTAD_results"
mkdir -p "${OUTPUT_DIR}"

# Extract unique resolutions
RESOLUTIONS=($(ls ${MATRIX_DIR}/*_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5 | xargs -n 1 basename | awk -F'_' '{print $3}' | sort -u))

# Select the current resolution based on the Slurm job array index
RESOLUTION=${RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}

echo "Processing resolution: ${RESOLUTION}"

# Match the treatment and control matrices for the current resolution
TARGET_MATRIX="${MATRIX_DIR}/O_LTHSC_${RESOLUTION}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"
CONTROL_MATRIX="${MATRIX_DIR}/Y_LTHSC_${RESOLUTION}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"

# Match the TAD domains file for the target condition
TAD_DOMAINS="${TAD_DIR}/O_LTHSC_${RESOLUTION}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_domains.bed"

# Check if files exist
if [[ -f "${TARGET_MATRIX}" && -f "${CONTROL_MATRIX}" && -f "${TAD_DOMAINS}" ]]; then
    echo "Target matrix: ${TARGET_MATRIX}"
    echo "Control matrix: ${CONTROL_MATRIX}"
    echo "TAD domains: ${TAD_DOMAINS}"

    # Run hicDifferentialTAD
    hicDifferentialTAD \
        --targetMatrix "${TARGET_MATRIX}" \
        --controlMatrix "${CONTROL_MATRIX}" \
        --tadDomains "${TAD_DOMAINS}" \
        --outFileNamePrefix "${OUTPUT_DIR}/O_vs_Y_LTHSC_hicDifferentialTAD_${RESOLUTION}" \
        --pValue 0.01 \
        --mode all \
        --modeReject all \
        --threads 4

    echo "hicDifferentialTAD completed for resolution: ${RESOLUTION}"
else
    echo "Error: Missing input files for resolution ${RESOLUTION}. Skipping."
fi

conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicDifferentialTAD_replicates_merged.sh