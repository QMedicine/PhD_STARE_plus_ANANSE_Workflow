#!/bin/bash
#SBATCH --job-name=convert_sra_to_fastq
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/convert_sra/convert_sra_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/convert_sra/convert_sra_%A_%a.err
#SBATCH --array=0-27 # 28 tasks total
#SBATCH --time=04:00:00 # Adjust time if needed, conversion can take a while
#SBATCH --cpus-per-task=2 # Total CPUs = 2 * 28 = 56 (within limit)
#SBATCH --mem=16G       # Total Memory = 16G * 28 = 448G (within 450G limit)

# Exit script on any error
set -e
set -u
set -o pipefail

echo "Starting SRA to FASTQ conversion job $SLURM_JOB_ID, task $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Allocated CPUs per task: $SLURM_CPUS_PER_TASK"
# echo "Allocated Memory per task: ${SLURM_MEM_PER_TASK}M" # This variable is not always set by Slurm, commented out. The request is --mem=16G per task.

# Load the required module
module load sra-toolkit/3.1.0

# --- Configuration ---
# Final output directory for FASTQ files
output_base_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data/part_2_files_MPP_etc"
# Base directory for temporary files (ensure this path exists and you have write permissions)
base_temp_dir="/vast/scratch/users/stacpoole.q/Temp/Itokawa_et_al_2022_SRA_files"
# Path to the list of SRA accession IDs
sra_list_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_SRR_list_p2.txt"
# --- End Configuration ---

# --- Preparations ---
# Ensure the final output directory exists
mkdir -p "$output_base_dir"

# Create a unique temporary directory for this specific task within the base temp dir
# NOTE: Assumes the base_temp_dir already exists! Create it manually if needed:
# mkdir -p /vast/scratch/users/stacpoole.q/Temp/Itokawa_et_al_2022_SRA_files
task_temp_dir="${base_temp_dir}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$task_temp_dir"
echo "Using temporary directory for this task: $task_temp_dir"

# Read the list of SRA accession IDs into an array
mapfile -t sra_ids < "$sra_list_file"

# Check if the array index is valid
if [ $SLURM_ARRAY_TASK_ID -ge ${#sra_ids[@]} ]; then
  echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is out of bounds for SRA list file $sra_list_file (size ${#sra_ids[@]})."
  # Clean up temp dir before exiting
  rm -rf "$task_temp_dir"
  exit 1
fi

# Use the SLURM_ARRAY_TASK_ID to select the correct accession
sra_id=${sra_ids[$SLURM_ARRAY_TASK_ID]}

# Construct the path to the input SRA file within its subdirectory in the output directory
input_sra_dir="$output_base_dir/$sra_id"
input_sra_file="$input_sra_dir/$sra_id.sra" # Common pattern for prefetch with --output-directory

# Check if the input SRA directory and file exist
if [ ! -d "$input_sra_dir" ]; then
    echo "Error: Input directory $input_sra_dir not found for SRA ID $sra_id."
    # Attempt to find .sra file directly in output_base_dir as a fallback
    input_sra_file_fallback="$output_base_dir/$sra_id.sra"
    if [ -f "$input_sra_file_fallback" ]; then
        echo "Found SRA file directly in $output_base_dir. Adjusting input path."
        input_sra_file="$input_sra_file_fallback"
        input_sra_dir="$output_base_dir" # Adjust input dir for logging/context if needed
    else
        echo "Also checked $input_sra_file_fallback and it was not found."
        rm -rf "$task_temp_dir" # Clean up temp dir before exiting
        exit 1
    fi
elif [ ! -f "$input_sra_file" ]; then
    # Check for common alternative naming conventions
    input_sra_file_alt1="$input_sra_dir/$sra_id"
    input_sra_file_alt2="$input_sra_dir/sra/$sra_id.sra" # Another possibility
     if [ -f "$input_sra_file_alt1" ]; then
        echo "Found SRA file as $input_sra_file_alt1. Adjusting input path."
        input_sra_file="$input_sra_file_alt1"
     elif [ -f "$input_sra_file_alt2" ]; then
        echo "Found SRA file as $input_sra_file_alt2. Adjusting input path."
        input_sra_file="$input_sra_file_alt2"
     else
        echo "Error: Input SRA file $input_sra_file not found within directory $input_sra_dir for SRA ID $sra_id."
        echo "Also checked $input_sra_file_alt1 and $input_sra_file_alt2 and they were not found."
        ls -lhR "$input_sra_dir" # List contents recursively for debugging
        rm -rf "$task_temp_dir" # Clean up temp dir before exiting
        exit 1
     fi
fi

# --- Run fasterq-dump ---
echo "Processing SRA ID: $sra_id"
echo "Input SRA file: $input_sra_file"
echo "Output directory for FASTQ: $output_base_dir"
echo "Temporary directory for conversion: $task_temp_dir"

echo "Running fasterq-dump..."
fasterq-dump \
    --outdir "$output_base_dir" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --progress \
    --split-files \
    --temp "$task_temp_dir" \
    "$input_sra_file"

dump_exit_code=$?
echo "fasterq-dump completed for $sra_id with exit code $dump_exit_code."

# --- Cleanup ---
# Clean up temporary directory regardless of exit code, but report if files might remain on failure
echo "Cleaning up temporary directory: $task_temp_dir"
rm -rf "$task_temp_dir"

if [ $dump_exit_code -eq 0 ]; then
    echo "Listing created FASTQ files in $output_base_dir for $sra_id:"
    ls -lh "$output_base_dir/${sra_id}"*.fastq || echo "No FASTQ files found matching ${sra_id}*.fastq"

    # --- Optional: Delete SRA file/directory after successful conversion ---
    # Uncomment the following lines ONLY if you are sure you want to delete
    # the original SRA data immediately after conversion to save space.
    # echo "Optional: Removing original SRA data directory $input_sra_dir"
    # rm -rf "$input_sra_dir"
    # echo "SRA data removal complete."
    # -----------------------------------------------------------------------

    echo "Task $SLURM_ARRAY_TASK_ID ($sra_id) finished successfully."
else
    echo "fasterq-dump failed for $sra_id with exit code $dump_exit_code."
    # Exit with the same error code
    exit $dump_exit_code
fi

# --- End Script ---

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/download_and_convert_SRA_files_for_Itokawa_et_al_2022_RNAseq_p2_convertSRA_to_fastq.sh