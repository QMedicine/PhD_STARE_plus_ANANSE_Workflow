#!/bin/bash
#SBATCH --job-name=fastqc_array          # Job name
#SBATCH --array=1-4                     # Array range (adjust if you have a different number of samples)
#SBATCH --cpus-per-task=8                # Number of CPUs per task
#SBATCH --mem=5G                        # Memory per job
#SBATCH --time=02:00:00                  # Time limit hrs:min:sec
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/fastqc_logs/array_job_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/fastqc_logs/array_job_%A_%a.err

# Load the FastQC module
module load fastqc

# Define directories
fastq_files_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data"
#scratch_output_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_RNAseq"
qc_output_dir="$fastq_files_dir/qc_reports/pre_trim"
tsv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_sample_information.tsv"

# Create output directory if it doesn't exist
mkdir -p "$qc_output_dir"

# Read the header line
header=$(head -n 1 "$tsv_file")

# Debugging output
echo "Header line: $header"

# Use awk to get the column indices for Read1, Read2, and Sample_ID
read1_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);  # Remove surrounding quotes
        if($i=="Read1") print i
    }
}')
read2_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);
        if($i=="Read2") print i
    }
}')
sample_id_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);
        if($i=="Sample_ID") print i
    }
}')

# Debugging outputs for column indices
echo "Read1 column: $read1_col"
echo "Read2 column: $read2_col"
echo "Sample_ID column: $sample_id_col"

# Check if column indices were found
if [[ -z "$read1_col" || -z "$sample_id_col" ]]; then
    echo "Error: Could not find Read1 or Sample_ID columns in header."
    exit 1
fi

# Read the specific sample info for this task
sample_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$tsv_file")

# Debugging output
echo "Sample line: $sample_line"

# Use awk to extract the fields from the sample line
read -r Sample_ID Read1 Read2 <<< $(echo "$sample_line" | awk -F'\t' -v read1_col="$read1_col" -v read2_col="$read2_col" -v sample_id_col="$sample_id_col" '{
    gsub(/^"|"$/, "", $read1_col);  # Remove quotes from Read1
    gsub(/^"|"$/, "", $read2_col);  # Remove quotes from Read2
    gsub(/^"|"$/, "", $sample_id_col);  # Remove quotes from Sample_ID
    printf "%s %s %s", $sample_id_col, $read1_col, $read2_col
}')

# Debugging outputs to check variables
echo "Sample_ID: $Sample_ID"
echo "Read1: $Read1"
echo "Read2: $Read2"

# Construct full paths to the reads
read1_file="$fastq_files_dir/$Read1"
read2_file="$fastq_files_dir/$Read2"

# Debugging outputs to check if the paths and variables are set correctly
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
    echo "Running FastQC on single-end data for sample '$Sample_ID'"
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read1_file"
else
    # Paired-end - check if Read2 file exists
    if [[ ! -f "$read2_file" ]]; then
        echo "Error: Read2 file '$read2_file' does not exist."
        exit 1
    fi
    echo "Running FastQC on paired-end data for sample '$Sample_ID'"
    # Run FastQC on both Read1 and Read2
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read1_file"
    fastqc -o "$qc_output_dir" -f fastq --noextract --threads 8 "$read2_file"
fi

# save and run on slurm using:
# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_fastqc.sh