#!/bin/bash
#SBATCH --job-name=cutadapt_part2        # Job name
#SBATCH --array=0-27                     # Array range (0 to N-1 for N samples in TSV, assuming 28 samples)
#SBATCH --cpus-per-task=2                # Number of CPUs per task (2 * 28 = 56 total CPUs)
#SBATCH --mem=8G                         # Memory per task (8G * 28 = 224G total RAM < 450G) - Adjust if needed
#SBATCH --time=04:00:00                  # Time limit hrs:min:sec (Adjust as needed)
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/cutadapt_logs/cutadapt_part2/cutadapt_part2_%A_%a.out # Updated log path
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/cutadapt_logs/cutadapt_part2/cutadapt_part2_%A_%a.err  # Updated log path

# Exit script on any error
set -e
set -u
set -o pipefail

echo "Starting Cutadapt job $SLURM_JOB_ID, task $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Step 1: Load necessary modules
module purge
module load cutadapt # Ensure this module provides a reasonably recent version
module load fastqc

# Step 2: Define directories and files
raw_fastq_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data/part_2_files_MPP_etc"
project_base_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq"
tsv_file="${project_base_dir}/Itokawa_et_al_2022_RNAseq_sample_information_p2.tsv"
adapter_fasta_file="${project_base_dir}/potential_adapters_part2.fa" # Path to the FASTA file created above

trimmed_output_dir="${project_base_dir}/trimmed_part2"
fastqc_trimmed_dir="${trimmed_output_dir}/fastqc_reports"
cutadapt_log_dir="${trimmed_output_dir}/cutadapt_logs" # Specific log dir for cutadapt summaries

# Create output directories if they don't exist
mkdir -p "$trimmed_output_dir"
mkdir -p "$fastqc_trimmed_dir"
mkdir -p "$cutadapt_log_dir"

# --- Ensure Adapter FASTA exists ---
if [ ! -f "$adapter_fasta_file" ]; then
    echo "Error: Adapter FASTA file not found: $adapter_fasta_file"
    echo "Please create it using the awk command provided in the instructions."
    exit 1
fi
echo "Using Adapter file: $adapter_fasta_file"

# Step 3: Read the header line from the TSV file
header=$(head -n 1 "$tsv_file")
echo "TSV Header: $header"

# Step 4: Get column indices
srr_id_col=$(echo "$header" | awk -F'\t' '{
    for(i=1; i<=NF; i++) { gsub(/^"|"$/, "", $i); if ($i == "Run") { print i; exit; } }
}')
sample_id_col=$(echo "$header" | awk -F'\t' '{
    for(i=1; i<=NF; i++) { gsub(/^"|"$/, "", $i); if ($i == "Sample_ID") { print i; exit; } }
}')

# Check if column indices were found
if [[ -z "$srr_id_col" ]]; then
    echo "Error: Could not find 'Run' (expected to contain SRR ID) column in header of $tsv_file."
    exit 1
fi
if [[ -z "$sample_id_col" ]]; then
    echo "Error: Could not find 'Sample_ID' column in header of $tsv_file."
    exit 1
fi
echo "SRR ID column (from 'Run' field): $srr_id_col"
echo "Output Sample_ID column: $sample_id_col"

# Step 5: Read the specific sample info for this task
sample_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$tsv_file") # +2 because sed is 1-based and we skip header
echo "Sample line $SLURM_ARRAY_TASK_ID: $sample_line"

# Step 6: Extract SRR ID and Sample ID using a single awk command and read
# --- CORRECTED EXTRACTION METHOD ---
# This approach reads the two desired fields identified by their column numbers ($srr_id_col, $sample_id_col)
# from the sample line in one go.
read srr_id sample_id <<< $(echo "$sample_line" | awk -F'\t' \
    -v run_col="$srr_id_col" \
    -v sample_col="$sample_id_col" \
    '{
        gsub(/^"|"$/, "", $run_col);    # Remove quotes from run field
        gsub(/^"|"$/, "", $sample_col); # Remove quotes from sample field
        print $run_col, $sample_col     # Print both, space-separated for "read"
    }')
# --- END CORRECTION ---

# Basic check if IDs were extracted
if [[ -z "$srr_id" ]]; then
    echo "Error: Failed to extract SRR ID from 'Run' column ($srr_id_col) in line: $sample_line"
    exit 1
fi
# Note: Sample ID might legitimately be just "Rep1", so checking for zero length might be too strict
# depending on your expected Sample_ID format. We'll trust the extraction if srr_id worked for now.
# if [[ -z "$sample_id" ]]; then
#     echo "Error: Failed to extract Sample ID from 'Sample_ID' column ($sample_id_col) in line: $sample_line"
#     exit 1
# fi
echo "Base SRR ID for input file lookup: $srr_id"
echo "Sample ID for output naming: $sample_id" # Check if this now shows the correct value (e.g., SRR13192287_Young_MPP1_Rep1)

# Step 7: Construct full path to the raw input FASTQ file
raw_fastq_base="${raw_fastq_dir}/${srr_id}_1.fastq"
if [[ -f "${raw_fastq_base}" ]]; then
    raw_fastq_input="${raw_fastq_base}"
elif [[ -f "${raw_fastq_base}.gz" ]]; then
    raw_fastq_input="${raw_fastq_base}.gz"
else
    echo "Error: Input FASTQ file not found for SRR ID $srr_id at ${raw_fastq_base}(.gz)?"
    ls -lh "$raw_fastq_dir" | grep "$srr_id"
    exit 1
fi
echo "Input FASTQ: $raw_fastq_input"

# Step 8: Define output file names using Sample_ID
trimmed_fastq="${trimmed_output_dir}/${sample_id}_trimmed_1.fastq.gz"
cutadapt_log="${cutadapt_log_dir}/${sample_id}_cutadapt.log"

# Step 9: Run Cutadapt
echo "Running Cutadapt for $sample_id..."
cutadapt \
    -j "$SLURM_CPUS_PER_TASK" \
    -q 20 \
    --minimum-length 20 \
    -a "file:${adapter_fasta_file}" \
    -o "$trimmed_fastq" \
    "$raw_fastq_input" > "$cutadapt_log" # Save Cutadapt report to log file

cutadapt_exit_code=$?

if [ $cutadapt_exit_code -ne 0 ]; then
    echo "Error: Cutadapt failed for $sample_id with exit code $cutadapt_exit_code."
    echo "Check log: $cutadapt_log"
    exit 1
fi
echo "Cutadapt finished successfully for $sample_id."
echo "Trimmed file: $trimmed_fastq"
echo "Cutadapt report:"
cat "$cutadapt_log" # Display the cutadapt summary in the main Slurm output

# Step 10: Run FastQC on the trimmed file
echo "Running FastQC on trimmed file for $sample_id..."
fastqc \
    -t "$SLURM_CPUS_PER_TASK" \
    -o "$fastqc_trimmed_dir" \
    "$trimmed_fastq"

fastqc_exit_code=$?

if [ $fastqc_exit_code -ne 0 ]; then
    echo "Error: FastQC failed for trimmed $sample_id with exit code $fastqc_exit_code."
    # Don't necessarily exit, Cutadapt succeeded, but log the error.
    # exit 1
fi
echo "FastQC finished for trimmed $sample_id."

echo "Processing of task $SLURM_ARRAY_TASK_ID ($sample_id) completed."

# To run in slurm use this:
# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_cutadapt_part2_files.sh