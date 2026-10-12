#!/bin/bash
#SBATCH --job-name=convertToH5
#SBATCH --array=0-35
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicConvertFormat_logs/hicConvertFormat_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicConvertFormat_logs/hicConvertFormat_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# Define directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files"
mkdir -p ${OUTPUT_DIR}

# Get list of .cool files
FILES=($(ls ${INPUT_DIR}/*.cool))

# Select the file for this array index
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename ${FILE} .cool)

# Convert .cool to .h5
hicConvertFormat \
    -m ${FILE} \
    --inputFormat cool \
    --outputFormat h5 \
    -o ${OUTPUT_DIR}/${BASENAME}.h5

# Deactivate the conda environment
conda deactivate

#sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicConvertFormat_cool_2_h5.sh