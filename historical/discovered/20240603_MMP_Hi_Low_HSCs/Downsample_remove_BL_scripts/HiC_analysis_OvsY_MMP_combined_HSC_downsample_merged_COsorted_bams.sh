#!/bin/bash
#SBATCH --job-name=downsampl_merged_COsorted_bams
#SBATCH --array=0-5                             # Adjust based on the number of samples
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampling_merged_log_files/Downsampling_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampling_merged_log_files/Downsampling_%A_%a.err

set -e  # Exit immediately if a command exits with a non-zero status

# Load required modules
module purge
module load samtools

# Define directories
BAM_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files"
OUTPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/downsampling_metrics"
DOWNSAMPLED_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/COSorted_downsampled_merged_bam_files"

# Create necessary directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DOWNSAMPLED_DIR}"

# Path to unique read counts file
UNIQUE_READ_COUNTS="${OUTPUT_DIR}/unique_read_counts.txt"

# Ensure unique_read_counts.txt exists
if [ ! -f "${UNIQUE_READ_COUNTS}" ]; then
    echo "Error: Unique read counts file not found at ${UNIQUE_READ_COUNTS}. Please generate it before running this script."
    exit 1
fi

# Get list of BAM files
BAM_FILES=("${BAM_DIR}"/*.bam)
TOTAL_SAMPLES=${#BAM_FILES[@]}

echo "Total number of BAM files found: ${TOTAL_SAMPLES}"

# Validate SLURM_ARRAY_TASK_ID
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${TOTAL_SAMPLES}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds number of samples (${TOTAL_SAMPLES}). Exiting."
    exit 1
fi

# Select the current BAM file based on the job array index
DEDUP_BAM="${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}"
SAMPLE=$(basename "${DEDUP_BAM}" .bam)

echo "Processing Sample: ${SAMPLE}"
echo "BAM File: ${DEDUP_BAM}"

# Find the minimum unique read count
MIN_UNIQUE_READS=$(awk 'NR==1 || $1 < min {min=$1} END {print min}' "${UNIQUE_READ_COUNTS}")
MIN_SAMPLE=$(awk -v min="${MIN_UNIQUE_READS}" '$1 == min {print $2}' "${UNIQUE_READ_COUNTS}")

echo "Minimum Unique Read Pairs: ${MIN_UNIQUE_READS} (Sample: ${MIN_SAMPLE})"

# Get unique read pairs for the current sample
CURRENT_SAMPLE_UNIQUE_READS=$(awk -v sample="${SAMPLE}" '$2 == sample {print $1}' "${UNIQUE_READ_COUNTS}")

# Validate unique read count
if [ -z "${CURRENT_SAMPLE_UNIQUE_READS}" ]; then
    echo "Error: Unique read count for sample ${SAMPLE} not found in ${UNIQUE_READ_COUNTS}. Exiting."
    exit 1
fi

echo "Unique Read Pairs for ${SAMPLE}: ${CURRENT_SAMPLE_UNIQUE_READS}"

# Determine if downsampling is needed
if [ "${CURRENT_SAMPLE_UNIQUE_READS}" -le "${MIN_UNIQUE_READS}" ]; then
    echo "Sample ${SAMPLE} has the minimum unique reads or fewer. No downsampling needed."
    # Copy the original BAM and index to downsampled directory
    cp "${DEDUP_BAM}" "${DOWNSAMPLED_DIR}/${SAMPLE}.bam"
    cp "${DEDUP_BAM}.bai" "${DOWNSAMPLED_DIR}/${SAMPLE}.bam.bai"
    
    # Generate a report indicating no downsampling
    REPORT_FILE="${OUTPUT_DIR}/${SAMPLE}_downsampling_report.txt"
    {
        echo "Sample: ${SAMPLE}"
        echo "Original Unique Read Pairs: ${CURRENT_SAMPLE_UNIQUE_READS}"
        echo "Downsampled Unique Read Pairs: ${CURRENT_SAMPLE_UNIQUE_READS}"
        echo "Downsampling Fraction: 1.000000"
        echo "Random Seed Used: N/A (No downsampling)"
    } > "${REPORT_FILE}"
    
    echo "Copied original BAM to downsampled directory and generated report at ${REPORT_FILE}"
    exit 0
fi

# Calculate downsampling fraction
DOWNSAMPLE_FRACTION=$(echo "scale=6; ${MIN_UNIQUE_READS} / ${CURRENT_SAMPLE_UNIQUE_READS}" | bc)

# Ensure DOWNSAMPLE_FRACTION is correctly formatted with leading zero
DOWNSAMPLE_FRACTION=$(printf "%.6f" "${DOWNSAMPLE_FRACTION}")

echo "Downsampling Fraction for ${SAMPLE}: ${DOWNSAMPLE_FRACTION}"

# Define the downsampled BAM file path
DOWNSAMPLED_BAM="${DOWNSAMPLED_DIR}/${SAMPLE}.downsampled.bam"

# Set a fixed random seed for reproducibility
RANDOM_SEED=43

# Downsample the BAM file using --subsample and --subsample-seed
echo "Downsampling ${SAMPLE} with subsample fraction: ${DOWNSAMPLE_FRACTION} and random seed: ${RANDOM_SEED}..."
samtools view -@ "${SLURM_CPUS_PER_TASK}" --subsample "${DOWNSAMPLE_FRACTION}" --subsample-seed "${RANDOM_SEED}" -b "${DEDUP_BAM}" -o "${DOWNSAMPLED_BAM}"

echo "Downsampling completed for ${SAMPLE}. Downsampled BAM: ${DOWNSAMPLED_BAM}"

# Index the downsampled BAM file
echo "Indexing downsampled BAM file for ${SAMPLE}..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${DOWNSAMPLED_BAM}"
echo "Indexing completed for ${DOWNSAMPLED_BAM}"

# Get downsampled unique read pairs
DOWNSAMPLED_READ_PAIRS=$(samtools view -c -F 1796 "${DOWNSAMPLED_BAM}")

# Generate downsampling report
REPORT_FILE="${OUTPUT_DIR}/${SAMPLE}_downsampling_report.txt"
{
    echo "Sample: ${SAMPLE}"
    echo "Original Unique Read Pairs: ${CURRENT_SAMPLE_UNIQUE_READS}"
    echo "Downsampled Unique Read Pairs: ${DOWNSAMPLED_READ_PAIRS}"
    echo "Downsampling Fraction: ${DOWNSAMPLE_FRACTION}"
    echo "Random Seed Used: ${RANDOM_SEED}"
} > "${REPORT_FILE}"

echo "Downsampling completed for ${SAMPLE}. Report saved to ${REPORT_FILE}."


# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_combined_HSC_downsample_merged_COsorted_bams.sh


