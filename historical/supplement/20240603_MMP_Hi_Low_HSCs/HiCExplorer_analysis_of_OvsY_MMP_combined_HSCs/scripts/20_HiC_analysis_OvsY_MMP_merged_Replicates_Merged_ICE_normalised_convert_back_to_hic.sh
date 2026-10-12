#!/bin/bash
#SBATCH --job-name=convert_2_hic
#SBATCH --array=0-5
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer/logs/hicConvertFormat_to_hic_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer/logs/hicConvertFormat_to_hic_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# This script is to convert the Hi-C h5 files for the MMP-merged, replicate-merged files that have been ICE normalized back to cool files in order to use the cool files for plotting heatmaps

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda init --user
conda activate hicexplorer

# Define directories
INPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer"
OUTPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer"
mkdir -p ${OUTPUT_DIR}

# Get list of .cool files
FILES=($(ls ${INPUT_DIR}/*.h5))

# Select the file for this array index
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename ${FILE} .h5)

# Convert .h5 to .hic
hicConvertFormat \
    -m ${FILE} \
    --inputFormat h5 \
    --outputFormat cool \
    -o ${OUTPUT_DIR}/${BASENAME}.cool

# Deactivate the conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/scripts/20_HiC_analysis_OvsY_MMP_merged_Replicates_Merged_ICE_normalised_convert_back_to_hic.sh
