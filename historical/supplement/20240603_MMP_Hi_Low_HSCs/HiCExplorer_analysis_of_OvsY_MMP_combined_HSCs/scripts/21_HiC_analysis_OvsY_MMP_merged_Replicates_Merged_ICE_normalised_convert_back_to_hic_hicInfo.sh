#!/bin/bash
#SBATCH --job-name=hicInfoJob       # Job name
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer/logs/21_hicInfo_%j.log      # Output log file
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer/logs/21_hicInfo_%j.err       # Error log file
#SBATCH --time=02:00:00              # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Load Miniconda3 module
module load miniconda3

# Activate the hicexplorer environment
conda activate hicexplorer

# Directories containing .cool files of the MMP merged replicate merged ICE normalized cool files
POSTCORR_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer"

# Create the output directory if it doesn't exist
OUTPUT_DIR="${POSTCORR_DIR}/hicInfo_results"
mkdir -p "$OUTPUT_DIR"

# Define the final output filename
final_output_file="${OUTPUT_DIR}/hicInfo_results.txt"

# Create or clear the final output file
> "$final_output_file"

# Create an array of pre-normalized and post-corrected files
post_files=("$POSTCORR_DIR"/*.cool)

# Loop directly through the file paths
for post_matrix in "${post_files[@]}"; do

    # Extract the base filename without extension for output naming
    base_filename=$(basename "$post_matrix" .cool)

    {
        # Append post-corrected information to the final output file
        # Use $base_filename here if you want cleaner logs
        echo "=== hicInfo for the ICE normalized cool file: $base_filename ==="
        hicInfo -m "$post_matrix"
        echo ""
        echo "--------------------------------------------"
    } >> "$final_output_file"
done

# Deactivate the conda environment after completion
conda deactivate