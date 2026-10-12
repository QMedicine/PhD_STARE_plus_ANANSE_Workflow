#!/bin/bash
#SBATCH --job-name=filter_bam_blacklist
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/logs/01_filter_bam_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/logs/01_filter_bam_%A_%a.err
#SBATCH --array=0-12
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# Define directories
BAM_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/filtered_bam_files"
TEMP_DIR="/vast/scratch/users/stacpoole.q/Temp_bam_files"
BLACKLIST="/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/filtered_bam_files_BL_removed"

# Ensure temporary directory exists
mkdir -p $TEMP_DIR

# Load necessary modules
module load samtools
module load bedtools

# Get the list of BAM files and select the current one
BAM_FILES=($BAM_DIR/*.bam)
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the base name of the BAM file
BASENAME=$(basename "$BAM_FILE" .bam)

# Step 1: Coordinate sort the BAM file
SORTED_BAM="$TEMP_DIR/${BASENAME}_sorted.bam"
samtools sort -@ 4 -o $SORTED_BAM $BAM_FILE

# Step 2: Filter out reads overlapping the blacklist regions
FILTERED_BAM="$TEMP_DIR/${BASENAME}_filtered.bam"
bedtools intersect -v -abam $SORTED_BAM -b $BLACKLIST > $FILTERED_BAM

# Step 3: Name sort the filtered BAM file and save the final output
FINAL_BAM="$OUTPUT_DIR/${BASENAME}_BL_removed.bam"
samtools sort -n -@ 4 -o $FINAL_BAM $FILTERED_BAM

# Clean up temporary files
rm -f $SORTED_BAM $FILTERED_BAM

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/filtered_bam_files/filter_bam_blacklist.sh