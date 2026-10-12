#!/bin/bash
#SBATCH --job-name=compare_raw_corrected
#SBATCH --array=0-35
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicCompareLogs/compare_raw_corrected_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicCompareLogs/compare_raw_corrected_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G

module load miniconda3
conda activate hicexplorer

# Directories for raw and corrected matrices
RAW_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files"
CORRECTED_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files/raw_v_corrected_comparison_plots"
mkdir -p ${OUTPUT_DIR}

# Get the list of raw and corrected files
RAW_FILES=($(ls ${RAW_DIR}/*.h5 | sort))
CORRECTED_FILES=($(ls ${CORRECTED_DIR}/*.h5 | sort))

# Select the current pair of files based on the Slurm array task ID
RAW_FILE=${RAW_FILES[$SLURM_ARRAY_TASK_ID]}
CORRECTED_FILE=${CORRECTED_FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${RAW_FILE}" .h5)

# Diagnostic plots for raw and corrected
hicCorrectMatrix diagnostic_plot \
    --matrix ${RAW_FILE} \
    --plotName ${OUTPUT_DIR}/${BASENAME}_raw_diagnostic.png

hicCorrectMatrix diagnostic_plot \
    --matrix ${CORRECTED_FILE} \
    --plotName ${OUTPUT_DIR}/${BASENAME}_corrected_diagnostic.png

# Heatmap visualization for raw and corrected with titles
hicPlotMatrix \
    --matrix ${RAW_FILE} \
    --title "${BASENAME} chr1:1-10265138 - Raw" \
    --region chr1:1-10265138 \
    --outFileName ${OUTPUT_DIR}/${BASENAME}_raw_heatmap_matrix.png

hicPlotMatrix \
    --matrix ${CORRECTED_FILE} \
    --title "${BASENAME} chr1:1-10265138 - Corrected" \
    --region chr1:1-10265138 \
    --outFileName ${OUTPUT_DIR}/${BASENAME}_corrected_heatmap_matrix.png

conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_raw_v_corrected_plots.sh