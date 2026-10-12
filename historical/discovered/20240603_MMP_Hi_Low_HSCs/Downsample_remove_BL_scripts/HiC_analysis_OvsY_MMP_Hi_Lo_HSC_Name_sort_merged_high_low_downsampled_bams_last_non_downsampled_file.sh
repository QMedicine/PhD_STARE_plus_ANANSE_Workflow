#!/bin/bash
#SBATCH --job-name=name_sort_downsampled_bams                  # Job name
#SBATCH --cpus-per-task=4                          # Number of CPUs per task
#SBATCH --mem=16G                                  # Memory per job
#SBATCH --time=2:00:00                             # Time limit (HH:MM:SS)
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/name_sort_downsampled_bam_log/NameSortBam_last_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/name_sort_downsampled_bam_log/NameSortBam_last_%A_%a.err

set -e  # Exit immediately if a command exits with a non-zero status

### **1. Load Required Modules**

# Clear any loaded modules and load samtools
module purge
module load samtools

### **2. Define Directories**

# Input directory containing the merged coordinate-sorted BAM files
INPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/COSorted_downsampled_merged_bam_files"

# Output directory for the name-sorted BAM files
OUTPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/Name_sorted_downsampled_merged_bam_files"

COORDINATE_SORTED_BAM="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/COSorted_downsampled_merged_bam_files/Y_LTHSC_MMP_combined_Rep1_COsorted.bam"
NAME_SORTED_BAM="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_bam_files/Name_sorted_downsampled_merged_bam_files/Y_LTHSC_MMP_combined_Rep1_Namesorted.downsampled.bam"


# Create output and log directories if they don't exist
mkdir -p "${OUTPUT_DIR}"

### **3. Perform SAMTOOLS Name Sort on last non-downsampled file**

echo "Name sorting last non-downsampled file ..."
samtools sort -n -@ "${SLURM_CPUS_PER_TASK}" -o "${NAME_SORTED_BAM}" "${COORDINATE_SORTED_BAM}"
echo "Name sorting completed


# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_Hi_Lo_HSC_Name_sort_merged_high_low_downsampled_bams_last_non_downsampled_file.sh
