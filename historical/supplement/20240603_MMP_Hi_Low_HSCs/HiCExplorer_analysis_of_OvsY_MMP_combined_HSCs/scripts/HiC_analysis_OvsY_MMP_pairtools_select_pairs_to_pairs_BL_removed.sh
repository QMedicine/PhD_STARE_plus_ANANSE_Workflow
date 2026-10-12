#!/bin/bash
#SBATCH --job-name=remove_blacklist
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/pairtools_remove_blacklist_logs/remove_blacklist_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/pairtools_remove_blacklist_logs/remove_blacklist_%A_%a.err
#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Load Miniconda and activate the environment
module load miniconda3
conda activate pairtools_py310

# Define directories and files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files"
OUTPUT_DIR="${INPUT_DIR}/Pairs_files_BL_removed"
CONDITIONS_FILE="/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges_as_pairtools_select_conditions.txt"
SCRIPT="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_pairtools_select_pairs_to_pairs_BL_removed.py"

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}/logs"

# Get all input files
FILES=(${INPUT_DIR}/*.nodups.pairs.filt.gz)

# Get the file corresponding to this job array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Run the Python script
python "${SCRIPT}" "${FILE}" "${OUTPUT_DIR}" "${CONDITIONS_FILE}"

# Deactivate the environment
conda deactivate

echo "Finished processing ${FILE}"

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_pairtools_select_pairs_to_pairs_BL_removed.sh