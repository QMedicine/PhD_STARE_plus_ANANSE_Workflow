#!/bin/bash
#SBATCH --job-name=merge_high_low_bams                # Job name
#SBATCH --array=0-5                               # Array indices for 6 samples (0-5)
#SBATCH --cpus-per-task=4                         # Number of CPUs per task
#SBATCH --mem=16G                                  # Memory per job
#SBATCH --time=4:00:00                             # Time limit (HH:MM:SS)
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/combine_bam_log/CombineBam_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/combine_bam_log/CombineBam_%A_%a.err

set -e  # Exit immediately if a command exits with a non-zero status

### **1. Load Required Modules**

# Clear any loaded modules and load samtools
module purge
module load samtools

### **2. Define Directories**

# Input directory containing the original COsorted BAM files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files"

# Output directory for the combined BAM files
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files"

# Create output and log directories if they don't exist
mkdir -p "${OUTPUT_DIR}"

### **3. Define Arrays for Cell Types and Replicates**

# Arrays mapping task IDs to cell types and replicates
CELL_TYPES=("O_LTHSC" "O_LTHSC" "O_LTHSC" "Y_LTHSC" "Y_LTHSC" "Y_LTHSC")
REPLICATES=("Rep1" "Rep2" "Rep3" "Rep1" "Rep2" "Rep3")

# Total number of samples
TOTAL_SAMPLES=${#CELL_TYPES[@]}

### **4. Validate SLURM_ARRAY_TASK_ID**

# Ensure the task ID is within the valid range
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${TOTAL_SAMPLES}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) exceeds the number of samples (${TOTAL_SAMPLES}). Exiting."
    exit 1
fi

### **5. Select Current Cell Type and Replicate Based on Task ID**

# Get the current cell type and replicate
CELL_TYPE=${CELL_TYPES[${SLURM_ARRAY_TASK_ID}]}
REPLICATE=${REPLICATES[${SLURM_ARRAY_TASK_ID}]}

echo "----------------------------------------"
echo "Processing Sample: ${CELL_TYPE}, Replicate: ${REPLICATE}"
echo "----------------------------------------"

### **6. Define Input and Output BAM File Paths**

# Define the input High and Low BAM file paths
HIGH_BAM="${INPUT_DIR}/${CELL_TYPE}_MMP_High_${REPLICATE}_COsorted.bam"
LOW_BAM="${INPUT_DIR}/${CELL_TYPE}_MMP_Low_${REPLICATE}_COsorted.bam"

# Check if input BAM files exist
if [ ! -f "${HIGH_BAM}" ]; then
    echo "Error: High BAM file not found at ${HIGH_BAM}. Exiting."
    exit 1
fi

if [ ! -f "${LOW_BAM}" ]; then
    echo "Error: Low BAM file not found at ${LOW_BAM}. Exiting."
    exit 1
fi

# Define the output combined BAM file path
COMBINED_BAM="${OUTPUT_DIR}/${CELL_TYPE}_MMP_combined_${REPLICATE}_COsorted.bam"

# Define the combined BAM index file path
COMBINED_BAI="${COMBINED_BAM}.bai"

### **7. Check if Combined BAM Already Exists**

# If the combined BAM already exists, skip merging
if [ -f "${COMBINED_BAM}" ]; then
    echo "Combined BAM file ${COMBINED_BAM} already exists. Skipping merging."
    exit 0
fi

### **8. Perform SAMTOOLS Merge**

echo "Merging ${HIGH_BAM} and ${LOW_BAM} into ${COMBINED_BAM}..."
samtools merge -@ "${SLURM_CPUS_PER_TASK}" -f "${COMBINED_BAM}" "${HIGH_BAM}" "${LOW_BAM}"
echo "Merging completed for ${CELL_TYPE}, Replicate: ${REPLICATE}"

### **9. Index the Combined BAM File**

echo "Indexing ${COMBINED_BAM}..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${COMBINED_BAM}"
echo "Indexing completed for ${COMBINED_BAM}"

### **10. Generate a Merge Report**

REPORT_FILE="${OUTPUT_DIR}/${CELL_TYPE}_MMP_combined_${REPLICATE}_merge_report.txt"
{
    echo "Sample: ${CELL_TYPE}"
    echo "Replicate: ${REPLICATE}"
    echo "Merged Files:"
    echo "${HIGH_BAM}"
    echo "${LOW_BAM}"
    echo "Combined BAM: ${COMBINED_BAM}"
} > "${REPORT_FILE}"

echo "Merge report saved to ${REPORT_FILE}"
echo "----------------------------------------"
echo "Processing of Sample: ${CELL_TYPE}, Replicate: ${REPLICATE} is complete."
echo "----------------------------------------"

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_Hi_Lo_HSC_merge_high_low_bams.sh