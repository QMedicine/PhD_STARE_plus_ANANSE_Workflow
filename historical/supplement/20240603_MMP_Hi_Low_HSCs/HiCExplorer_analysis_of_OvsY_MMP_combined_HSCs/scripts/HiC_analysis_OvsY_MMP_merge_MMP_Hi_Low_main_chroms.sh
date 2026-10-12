#!/bin/bash
#SBATCH --job-name=merge_MMP_cool
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/merge_MMP_cool_logs/merge_MMP_cool_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/merge_MMP_cool_logs/merge_MMP_cool_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G

# Directory containing the .cool files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/Cool_files_BL_removed_main_chroms_only"
OUTPUT_DIR="${INPUT_DIR}/OandY_LTHSC_MMP_combined_Cool_files_main_chroms_only"
mkdir -p ${OUTPUT_DIR}

# Load Miniconda
module load miniconda3
conda init --user
conda activate cooler-env

# Navigate to the directory with .cool files
cd ${INPUT_DIR}

# Check if there are .cool files in the directory
if ! ls *.cool 1> /dev/null 2>&1; then
    echo "No .cool files found in ${INPUT_DIR}. Exiting."
    conda deactivate
    exit 1
fi

# Temporary file to store groups
GROUP_FILE=$(mktemp)

echo "Parsing and grouping files..."

for FILE in *.cool; do
    BASENAME=$(basename "$FILE")
    
    # Extract components
    CONDITION=$(echo "$BASENAME" | awk -F'_' '{print $1 "_" $2}')
    REPLICATE=$(echo "$BASENAME" | awk -F'_' '{print $5}')
    MMP_STATUS=$(echo "$BASENAME" | awk -F'_' '{print $3 "_" $4}')
    RESOLUTION=$(echo "$BASENAME" | awk -F'_' '{print $NF}' | sed 's/\.cool//')
    
    # Save grouping information to a temporary file
    echo "${CONDITION} ${REPLICATE} ${RESOLUTION} ${FILE}" >> ${GROUP_FILE}
done

# Process the groups
echo "Grouping complete. Processing merges..."
sort ${GROUP_FILE} | while read CONDITION REPLICATE RESOLUTION FILE; do
    # Find matching files with the same CONDITION, REPLICATE, and RESOLUTION
    MATCHING_FILES=$(grep -E "^${CONDITION} ${REPLICATE} ${RESOLUTION} " ${GROUP_FILE} | awk '{print $4}')
    
    # Count matching files
    NUM_FILES=$(echo "${MATCHING_FILES}" | wc -l)
    echo "Group: ${CONDITION}_${REPLICATE}_${RESOLUTION}"
    echo "  Files: ${MATCHING_FILES}"
    
    # Ensure there are exactly two files to merge
    if [[ ${NUM_FILES} -ne 2 ]]; then
        echo "  Warning: Expected 2 files, but found ${NUM_FILES}. Skipping this group."
        continue
    fi

    # Define output file name
    OUTPUT_FILE="${OUTPUT_DIR}/${CONDITION}_${REPLICATE}_${RESOLUTION}_MMP_combined_BL_removed.cool"
    
    # Merge files using cooler merge
    echo "Merging files: ${MATCHING_FILES}"
    if ! cooler merge ${OUTPUT_FILE} ${MATCHING_FILES}; then
        echo "Error merging files for group ${CONDITION}_${REPLICATE}_${RESOLUTION}. Skipping."
        continue
    fi
    
    echo "Successfully merged files into ${OUTPUT_FILE}"
done

# Cleanup
rm -f ${GROUP_FILE}

echo "Merging complete. Merged files are in ${OUTPUT_DIR}."

# Deactivate Conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/scripts/HiC_analysis_OvsY_MMP_merge_MMP_Hi_Low_main_chroms.sh

