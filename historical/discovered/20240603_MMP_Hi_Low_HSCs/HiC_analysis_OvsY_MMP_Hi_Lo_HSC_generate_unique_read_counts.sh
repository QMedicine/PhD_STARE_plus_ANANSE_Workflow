#!/bin/bash

# Directories
BAM_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files"
OUTPUT_DIR="${BAM_DIR}/downsampling_metrics"
mkdir -p "${OUTPUT_DIR}"

UNIQUE_READ_COUNTS="${OUTPUT_DIR}/unique_read_counts.txt"
> "${UNIQUE_READ_COUNTS}"

# Load samtools module
module purge
module load samtools

# Loop over BAM files
for BAM in "${BAM_DIR}"/*.bam; do
    SAMPLE=$(basename "${BAM}" .bam)
    echo "Processing ${SAMPLE}..."

    # Calculate number of unique read pairs
    UNIQUE_READ_PAIRS=$(samtools view -c -F 1796 "${BAM}")
    echo -e "${UNIQUE_READ_PAIRS}\t${SAMPLE}" >> "${UNIQUE_READ_COUNTS}"
done

echo "Unique read counts have been saved to ${UNIQUE_READ_COUNTS}"

# chmod +x /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_Hi_Lo_HSC_generate_unique_read_counts.sh