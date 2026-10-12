#!/bin/bash

# Input directory containing the .diff_tad files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files/hicFindTADs_TADs/hicDifferentialTAD_results"

# Loop through all *_rejected.diff_tad files
for file in ${INPUT_DIR}/*_rejected.diff_tad; do
    # Extract the base name without extension
    base_name=$(basename "${file}" _rejected.diff_tad)

    # Define the output .bed file
    output_file="${INPUT_DIR}/${base_name}_rejected.bed"

    # Convert to BED format by selecting the first 6 columns
    awk '!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "${file}" > "${output_file}"

    # Print a message for each processed file
    echo "Converted ${file} to ${output_file}"
done

echo "All .diff_tad files have been converted to .bed format."

# Loop through all *_accepted.diff_tad files
for file in ${INPUT_DIR}/*_accepted.diff_tad; do
    # Extract the base name without extension
    base_name=$(basename "${file}" _accepted.diff_tad)

    # Define the output .bed file
    output_file="${INPUT_DIR}/${base_name}_accepted.bed"

    # Convert to BED format by selecting the first 6 columns
    awk '!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "${file}" > "${output_file}"

    # Print a message for each processed file
    echo "Converted ${file} to ${output_file}"
done

echo "All .diff_tad files have been converted to .bed format."
