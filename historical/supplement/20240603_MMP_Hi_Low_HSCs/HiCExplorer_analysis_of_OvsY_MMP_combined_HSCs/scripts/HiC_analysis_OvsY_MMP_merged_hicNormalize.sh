#!/bin/bash
#SBATCH --job-name=normalize_h5
#SBATCH --array=0-35
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicNormalize_logs/hicNormalize_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicNormalize_logs/hicNormalize_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# Input and output directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files"
mkdir -p ${OUTPUT_DIR}

# List of all .h5 files in the input directory
FILES=($(ls ${INPUT_DIR}/*.h5))

# Select the file for the current array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "${FILE}" .h5)

# Normalize the matrix using hicNormalize
hicNormalize \
    --matrices "${FILE}" \
    --normalize smallest \
    --outFileName "${OUTPUT_DIR}/${BASENAME}_hicNormalize.h5"

# Deactivate conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicNormalize.sh