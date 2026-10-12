#!/bin/bash
#SBATCH --job-name=fastqc_array          # Job name
#SBATCH --array=1-36                     # Array range (adjust if you have a different number of samples)
#SBATCH --cpus-per-task=8                # Number of CPUs per task
#SBATCH --mem=20G                        # Memory per job
#SBATCH --time=02:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/array_job_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/array_job_%A_%a.err

# Load the FastQC module
module load fastqc

# Define directories
fastq_files_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq/raw_data"
scratch_output_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq"
qc_output_dir="$scratch_output_dir/qc_reports/pre_trim"
csv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/sample_information.csv"

# Create output directory if it doesn't exist
mkdir -p "$qc_output_dir"

# Get the column indices for Read1 and Read2
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

# Run FastQC on single-end or paired-end based on Read2
if [[ -z "$Read2" || "$Read2" == "NA" ]]; then
    # Single-end
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read1_file"
else
    # Paired-end - check if Read2 file exists
    if [[ ! -f "$read2_file" ]]; then
        echo "Error: Read2 file '$read2_file' does not exist."
        exit 1
    fi
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read1_file"
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read2_file"
fi


# Save this script as /vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/Itokawa_et_al_2022_ATACseq_job_array_fastqc.sh and run from a slurm login through terminal
# Use ssh slurm-login to login and enter password
# Run job using sbatch Itokawa_et_al_2022_ATACseq_job_array.sh