#!/bin/bash
#SBATCH --job-name=ATAC_job_array
#SBATCH --cpus-per-task=32
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --array=1-36
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/array_job_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/array_job_%A_%a.err

# Load necessary modules
module load bowtie2/2.5.3
module load picard-tools/3.1.1
module load MACS3/3.0.1 
module load openjdk/21.0.2
module load samtools/1.22.1
module load bedtools/2.31.1

# Set PICARD_HOME environment variable
export PICARD_HOME=/stornext/System/data/software/rhel/9/base/bioinf/picard-tools/3.1.1

# Define directories
fastq_files_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq/raw_data"
scratch_output_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq"
qc_output_dir="$scratch_output_dir/qc_reports/pre_trim"
project_output_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq"
genome_dir="/vast/projects/Immune_Aging/Reference_genomes/mm39/Bowtie2"
csv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/sample_information.csv"

g_size=2654621783 #GRCm39 effective genome size

# Create the various output directories if they don't exist
mkdir -p "$scratch_output_dir/aligned"
mkdir -p "$scratch_output_dir/deduped"
mkdir -p "$scratch_output_dir/deduped/metrics"
mkdir -p "$project_output_dir/filtered"
mkdir -p "$project_output_dir/macs3_peaks"

# Get the column indices for Read1, Read2, Sample_ID
header=$(head -n 1 "$csv_file")
read1_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Read1" | cut -d: -f1)
read2_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Read2" | cut -d: -f1)
sample_id_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Sample_ID" | cut -d: -f1)

# Read the specific sample info for this task
sample_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$csv_file")  # +1 because sed is 1-indexed

# Extract Read1, Read2, and Sample_ID fields from the CSV line
Read1=$(echo "$sample_line" | awk -v col=$read1_col 'BEGIN {FS=","} {print $col}' | sed 's/^"\|"$//g')
Read2=$(echo "$sample_line" | awk -v col=$read2_col 'BEGIN {FS=","} {print $col}' | sed 's/^"\|"$//g')
Sample_ID=$(echo "$sample_line" | awk -v col=$sample_id_col 'BEGIN {FS=","} {print $col}' | sed 's/^"\|"$//g')

# Construct full path to the reads
read1_file="$fastq_files_dir/$Read1"
read2_file="$fastq_files_dir/$Read2"

# Debugging outputs to check if the paths and variables are set correctly
echo "Sample_ID: $Sample_ID"
echo "Read1 file: $read1_file"
echo "Read2 file: $read2_file"

# Check if the Read1 file exists
if [[ ! -f "$read1_file" ]]; then
    echo "Error: Read1 file '$read1_file' does not exist."
    exit 1
fi

# Construct file names for output files
aligned_bam="$scratch_output_dir/aligned/${Sample_ID}_aligned.bam"
sorted_bam="$scratch_output_dir/aligned/${Sample_ID}_sorted.bam"
dedup_bam="$scratch_output_dir/deduped/${Sample_ID}_dedup.bam"
dedup_metrics="$scratch_output_dir/deduped/metrics/${Sample_ID}_dedup_metrics.txt"
filtered_bam="$project_output_dir/filtered/${Sample_ID}_filtered.bam"
peak_dir="$project_output_dir/macs3_peaks"

# Run Bowtie2 and other commands for single-end or paired-end reads
if [[ -z "$Read2" || "$Read2" == "NA" ]]; then
    # Single-end
#    bowtie2 -p 32 --very-sensitive -x "$genome_dir/mm39" -U "$read1_file" | \
#    samtools view -bS - > "$aligned_bam" && \
#    samtools sort -o "$sorted_bam" "$aligned_bam" && \
#    java -jar $PICARD_HOME/picard.jar MarkDuplicates I="$sorted_bam" O="$dedup_bam" M="$dedup_metrics" REMOVE_DUPLICATES=true && \
    samtools view -@ 32 -h -q 10 -F 4 "$dedup_bam" | grep -v 'chrM' | samtools view -@ 32 -b -o "$filtered_bam" && \
    macs3 callpeak -t "$filtered_bam" -g $g_size --nomodel --shift -100 --extsize 200 --keep-dup all -p 0.01 -n "$Sample_ID" --outdir "$peak_dir"

else
    # Paired-end
#    bowtie2 -p 32 --very-sensitive -x "$genome_dir/mm39" -1 "$read1_file" -2 "$read2_file" | \
#    samtools view -bS - > "$aligned_bam" && \
#    samtools sort -o "$sorted_bam" "$aligned_bam" && \
#    java -jar $PICARD_HOME/picard.jar MarkDuplicates I="$sorted_bam" O="$dedup_bam" M="$dedup_metrics" REMOVE_DUPLICATES=true && \
    samtools view -@ 32 -h -q 10 -F 4 "$dedup_bam" | grep -v 'chrM' | samtools view -@ 32 -b -o "$filtered_bam" && \
    macs3 callpeak -t "$filtered_bam" -f BAMPE -g $g_size --keep-dup all -p 0.01 -n "$Sample_ID" --outdir "$peak_dir"
fi

echo "Processing of $Sample_ID completed."

# Save this script as /vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/Itokawa_et_al_2022_ATACseq_job_array.sh and run from a slurm login through terminal
# Use ssh slurm-login to login and enter password
# Run job using sbatch Itokawa_et_al_2022_ATACseq_job_array.sh