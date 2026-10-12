#!/bin/bash

#SBATCH --job-name=count_reads
#SBATCH --array=0-11  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=1
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/count_reads_log/count_reads_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/count_reads_log/count_reads_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=24G

# Load samtools module if necessary
module load samtools

INPUTDIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/filtered_bam_files_BL_removed"

# Create or clear the output file if it doesn't exist
output_file="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/filtered_bam_files_BL_removed/interactions_per_chromosome.txt"
if [ ! -f "$output_file" ]; then
    echo -e "BAM_File\tChromosome\tCount" > "$output_file"
fi

# Get the BAM file corresponding to the array task ID
bam_files=(${INPUTDIR}/*.bam)
bam_file="${bam_files[$SLURM_ARRAY_TASK_ID]}"

# Check if the BAM file exists
if [[ ! -f "$bam_file" ]]; then
    echo "BAM file $bam_file does not exist." >> "$output_file"
    exit 1
fi

echo "Processing $bam_file..."

# Count reads for chromosomes 1-19, X, and Y using awk
samtools view "$bam_file" | awk '{
    if ($3 ~ /^chr[1-19XY]$/) counts[$3]++
} END {
    for (chr in counts) {
        print "'$bam_file'" "\t" chr "\t" counts[chr]
    }
}' >> "$output_file"

echo "Completed processing for $bam_file."



# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/count_reads_per_chromosome.sh
