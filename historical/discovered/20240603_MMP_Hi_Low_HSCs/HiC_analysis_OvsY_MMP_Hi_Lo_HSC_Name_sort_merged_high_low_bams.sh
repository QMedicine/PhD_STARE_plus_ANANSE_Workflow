#!/bin/bash
#SBATCH --job-name=name_sort_bams                  # Job name
#SBATCH --array=0-5                                # Array indices for 6 samples (0-5)
#SBATCH --cpus-per-task=4                          # Number of CPUs per task
#SBATCH --mem=16G                                  # Memory per job
#SBATCH --time=2:00:00                             # Time limit (HH:MM:SS)
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/name_sort_bam_log/NameSortBam_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/name_sort_bam_log/NameSortBam_%A_%a.err

set -e  # Exit immediately if a command exits with a non-zero status

### **1. Load Required Modules**

# Clear any loaded modules and load samtools
module purge
module load samtools

### **2. Define Directories**

# Input directory containing the merged coordinate-sorted BAM files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files"

# Output directory for the name-sorted BAM files
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files/Name_sorted_merged_high_low_bam_files"

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

# Define the input coordinate-sorted BAM file path
COORDINATE_SORTED_BAM="${INPUT_DIR}/${CELL_TYPE}_MMP_combined_${REPLICATE}_COsorted.bam"

# Check if input BAM file exists
if [ ! -f "${COORDINATE_SORTED_BAM}" ]; then
    echo "Error: Coordinate-sorted BAM file not found at ${COORDINATE_SORTED_BAM}. Exiting."
    exit 1
fi

# Define the output name-sorted BAM file path
NAME_SORTED_BAM="${OUTPUT_DIR}/${CELL_TYPE}_MMP_combined_${REPLICATE}_Namesorted.bam"

# Define the name-sorted BAM index file path
NAME_SORTED_BAI="${NAME_SORTED_BAM}.bai"

### **7. Check if Name-Sorted BAM Already Exists**

# If the name-sorted BAM already exists, skip sorting
if [ -f "${NAME_SORTED_BAM}" ]; then
    echo "Name-sorted BAM file ${NAME_SORTED_BAM} already exists. Skipping sorting."
    exit 0
fi

### **8. Perform SAMTOOLS Name Sort**

echo "Name sorting ${COORDINATE_SORTED_BAM} into ${NAME_SORTED_BAM}..."
samtools sort -n -@ "${SLURM_CPUS_PER_TASK}" -o "${NAME_SORTED_BAM}" "${COORDINATE_SORTED_BAM}"
echo "Name sorting completed for ${CELL_TYPE}, Replicate: ${REPLICATE}"

### **9. Index the Name-Sorted BAM File (Optional)**

echo "Indexing ${NAME_SORTED_BAM}..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${NAME_SORTED_BAM}"
echo "Indexing completed for ${NAME_SORTED_BAM}"

### **10. Generate a Sort Report**

REPORT_FILE="${OUTPUT_DIR}/${CELL_TYPE}_MMP_combined_${REPLICATE}_name_sort_report.txt"
{
    echo "Sample: ${CELL_TYPE}"
    echo "Replicate: ${REPLICATE}"
    echo "Input BAM: ${COORDINATE_SORTED_BAM}"
    echo "Output BAM: ${NAME_SORTED_BAM}"
} > "${REPORT_FILE}"

echo "Sort report saved to ${REPORT_FILE}"
echo "----------------------------------------"
echo "Processing of Sample: ${CELL_TYPE}, Replicate: ${REPLICATE} is complete."
echo "----------------------------------------"

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_Hi_Lo_HSC_Name_sort_merged_high_low_bams.sh