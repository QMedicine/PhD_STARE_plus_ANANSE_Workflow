#!/bin/bash
#SBATCH --job-name=remove_BL_regions
#SBATCH --array=0-19
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/hicAdjustMatrix_logs_merged_reps/hicNormalize_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/hicAdjustMatrix_logs_merged_reps/hicNormalize_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# Input and output directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp/cool_files_1000bp_BL_removed"
mkdir -p ${OUTPUT_DIR}

# List of all .cool files in the input directory
FILES=($(ls ${INPUT_DIR}/*.cool))

# Select the file for the current array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${FILE}" .cool)
echo "Input file: ${FILE}"
echo "Basename: ${BASENAME}"

# Normalize the matrix using hicNormalize
hicAdjustMatrix \
    --matrix "${FILE}" \
    --regions /vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed \
    --action remove \
    --outFileName "${OUTPUT_DIR}/${BASENAME}_BL_removed.cool"

echo "File saved to:${OUTPUT_DIR}/${BASENAME}_BL_removed.cool"

# Deactivate conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_hicAdjustMatrix_cool_1000bp_to_BL_removed.sh