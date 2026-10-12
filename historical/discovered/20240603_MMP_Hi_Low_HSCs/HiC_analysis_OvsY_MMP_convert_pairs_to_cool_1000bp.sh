#!/bin/bash
#SBATCH --job-name=pairs_to_cool_1000bp
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/pairs_to_cool_1000bp_logs/convert_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/pairs_to_cool_1000bp_logs/convert_%A_%a.err
#SBATCH --array=0-19
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load Miniconda
module load miniconda3
conda activate cooler-env

# Define variables
PAIRS_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files"
CHROM_SIZES="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/mm39.chrom.sizes"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp"
RESOLUTIONS=("1000")

# Make output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Get the file to process based on the array index
FILES=(${PAIRS_DIR}/*.nodups.pairs.filt.gz)
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#FILES[@]} ]; then
    echo "Array index out of bounds. Exiting."
    exit 1
fi
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename ${FILE} .nodups.pairs.filt.gz)

# Loop through resolutions and generate cool files
for RES in "${RESOLUTIONS[@]}"; do
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_${RES}.cool"
    echo "Processing ${FILE} at resolution ${RES}..."
    if ! cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${CHROM_SIZES}:${RES} ${FILE} ${OUTPUT_FILE}; then
        echo "Failed to process ${FILE} at resolution ${RES}"
    else
        echo "Successfully processed ${FILE} at resolution ${RES}"
    fi
done

echo "Conversion for ${BASENAME} completed."

conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_convert_pairs_to_cool_1000bp.sh