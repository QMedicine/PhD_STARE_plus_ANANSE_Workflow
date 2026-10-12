#!/bin/bash

# Combine all PNG plots into a single PDF
module load ImageMagick/7.1.1

# Directory containing the plots
PLOT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_diag_plots"
OUTPUT_FILE="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_diag_plots/diagnostic_plots.pdf"

OUTPUT_ZIP="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_diag_plots/hicCorrectMatrix_labeled_plots.zip"

# Temporary directory for labeled images
LABEL_DIR="${PLOT_DIR}/labeled_plots"
mkdir -p ${LABEL_DIR}

# Add labels below each plot
for FILE in ${PLOT_DIR}/*.png; do
    BASENAME=$(basename "${FILE}" .png)  # Extract the base name without extension
    magick \
        "${FILE}" \
        -gravity south \
        -background white -splice 0x50 \
        -pointsize 14 \
        -font "FreeSans" \
        -annotate +0+10 "${BASENAME}" \
        "${LABEL_DIR}/${BASENAME}_labeled.png"
done

# Combine labeled plots into a single PDF
magick ${LABEL_DIR}/*.png ${OUTPUT_FILE}

# Compress labeled plots into a zip file
zip -r ${OUTPUT_ZIP} ${LABEL_DIR}

# Cleanup temporary labeled images
rm -r ${LABEL_DIR}
