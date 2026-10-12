#!/bin/bash
#SBATCH --job-name=hicInfoJob       # Job name
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/hicInfoJob_logs/hicInfo_%j.log      # Output log file
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/hicInfoJob_logs/hicInfo_%j.err       # Error log file
#SBATCH --time=02:00:00              # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Load Miniconda3 module
module load miniconda3

# Activate the hicexplorer environment
conda activate hicexplorer

# Directories containing .h5 files
PRENORM_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files"
POSTCORR_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files"

# Create the output directory if it doesn't exist
OUTPUT_DIR="${POSTCORR_DIR}/hicInfo_compare_pre_post_correction_results"
mkdir -p "$OUTPUT_DIR"

# Define the final output filename
final_output_file="${OUTPUT_DIR}/hicInfo_comparison_pre_vs_post_correction.txt"

# Create or clear the final output file
> "$final_output_file"

# Create an array of pre-normalized and post-corrected files
pre_files=("$PRENORM_DIR"/*.h5)
post_files=("$POSTCORR_DIR"/*.h5)

# Check if both directories have the same number of files
if [ ${#pre_files[@]} -ne ${#post_files[@]} ]; then
    echo "Error: The number of pre-normalized and post-corrected files do not match."
    exit 1
fi

# Loop through the files by index
for i in "${!pre_files[@]}"; do
    pre_matrix="${pre_files[$i]}"
    post_matrix="${post_files[$i]}"

    # Extract the base filename without extension for output naming
    base_filename=$(basename "$pre_matrix" .h5)

    # Append pre-normalized information to the final output file
    {
        echo "=== Pre-normalized: $pre_matrix ==="
        hicInfo -m "$pre_matrix"
        echo ""
        
        # Append post-corrected information to the final output file
        echo "=== Post-corrected: $post_matrix ==="
        hicInfo -m "$post_matrix"
        echo ""
        echo "--------------------------------------------"
    } >> "$final_output_file"
done

# Deactivate the conda environment after completion
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/scripts/HiC_analysis_OvsY_MMP_merged_replicates_merged_hicInfo_pre_post_corrected.sh
