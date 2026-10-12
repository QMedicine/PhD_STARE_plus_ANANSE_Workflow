#!/bin/bash
#SBATCH --job-name=fastqc_part2         # Job name
#SBATCH --array=0-27                    # Array range (0 to N-1 for N samples in the list file)
#SBATCH --cpus-per-task=2               # CPUs per task (2 * 28 = 56 total CPUs)
#SBATCH --mem=5G                        # Memory per task (5G * 28 = 140G total RAM < 450G)
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec (Adjust as needed)
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/fastqc_logs/fastqc_part2/fastqc_part2_%A_%a.out # Updated log path
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/fastqc_logs/fastqc_part2/fastqc_part2_%A_%a.err  # Updated log path

# Exit script on any error
set -e
set -u
set -o pipefail

echo "Starting FastQC job $SLURM_JOB_ID, task $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Load the FastQC module (adjust version if necessary)
module load fastqc

# --- Configuration ---
# Directory containing the FASTQ files
fastq_files_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data/part_2_files_MPP_etc"
# Directory to save FastQC reports
qc_output_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/fastqc_reports/part2_pre_trim" # New output subdir
# List of SRA IDs used for file naming
sra_list_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_SRR_list_p2.txt"
# --- End Configuration ---

# Create output directory if it doesn't exist
mkdir -p "$qc_output_dir"

# Read the list of SRA accession IDs into an array
mapfile -t sra_ids < "$sra_list_file"

# Check if the array index is valid
if [ $SLURM_ARRAY_TASK_ID -ge ${#sra_ids[@]} ]; then
  echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is out of bounds for SRA list file $sra_list_file (size ${#sra_ids[@]})."
  exit 1
fi

# Use the SLURM_ARRAY_TASK_ID to select the correct accession
sra_id=${sra_ids[$SLURM_ARRAY_TASK_ID]}

# Construct the expected FASTQ filename (assuming single-end with _1 suffix)
fastq_file="${fastq_files_dir}/${sra_id}_1.fastq"

echo "Processing SRA ID: $sra_id"
echo "Expecting FASTQ file: $fastq_file"
echo "Output directory for report: $qc_output_dir"

# Check if the FASTQ file exists
if [[ ! -f "$fastq_file" ]]; then
    # Check if it might be gzipped
    if [[ -f "${fastq_file}.gz" ]]; then
        fastq_file="${fastq_file}.gz"
        echo "Found gzipped file: $fastq_file"
    else
      echo "Error: FASTQ file '$fastq_file' (or .gz version) not found."
      ls -lh "$fastq_files_dir" | grep "$sra_id" # List matching files for debugging
      exit 1
    fi
fi

# Run FastQC
echo "Running FastQC on $fastq_file ..."
fastqc \
    -o "$qc_output_dir" \
    -f fastq \
    --noextract \
    --threads "$SLURM_CPUS_PER_TASK" \
    "$fastq_file"

fastqc_exit_code=$?

if [ $fastqc_exit_code -eq 0 ]; then
    echo "FastQC finished successfully for $sra_id."
else
    echo "FastQC failed for $sra_id with exit code $fastqc_exit_code."
    exit $fastqc_exit_code
fi

echo "Task $SLURM_ARRAY_TASK_ID ($sra_id) complete."


# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_fastqc_part2_files.sh