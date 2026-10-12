#!/bin/bash

# This script generates individual Slurm submission scripts for STAR alignment.

# Exit on error
set -e
set -u
set -o pipefail

# --- Configuration ---
TRIMMED_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed_part2"
STAR_INDEX="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/STAR_index_60bp"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams"
# Use the metrics dir for Slurm .out/.err logs
LOG_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_alignment_log/part2_MPP_files_etc"
# Directory to store the generated .sh scripts
SCRIPT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/star_submission_scripts_part2"

# --- Resource Requests for Each Individual Job ---
CPUS_PER_JOB=32
MEM_PER_JOB="40G" # Using 40GB as requested (corrected from MB)
TIME_LIMIT="48:00:00"
SORT_RAM_LIMIT=30000000000 # 30 Billion Bytes (~30GB)

# --- Script Logic ---
echo "Generating STAR submission scripts..."

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${SCRIPT_DIR}"

# Find all trimmed fastq files
# Using find and while read is safer for filenames than simple globbing
find "${TRIMMED_DIR}" -maxdepth 1 -name '*_trimmed_1.fastq.gz' | sort | while IFS= read -r FASTQ_FILE; do

    # Extract the sample ID (part before _trimmed_1.fastq.gz)
    SAMPLE_ID=$(basename "${FASTQ_FILE}" "_trimmed_1.fastq.gz")

    # Define the name for the individual slurm script
    SLURM_SCRIPT_NAME="${SCRIPT_DIR}/submit_star_${SAMPLE_ID}.sh"

    echo "Creating script for Sample: ${SAMPLE_ID} -> ${SLURM_SCRIPT_NAME}"

    # Create the individual Slurm script using a "here document"
    cat << EOF > "${SLURM_SCRIPT_NAME}"
#!/bin/bash
#SBATCH --job-name=STAR_${SAMPLE_ID}
#SBATCH --output=${LOG_DIR}/STAR_${SAMPLE_ID}.out
#SBATCH --error=${LOG_DIR}/STAR_${SAMPLE_ID}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_JOB}
#SBATCH --mem=${MEM_PER_JOB}
#SBATCH --time=${TIME_LIMIT}

# Exit on error
set -e
set -u
set -o pipefail

echo "Starting STAR alignment for sample: ${SAMPLE_ID}"
echo "Running on node: \$(hostname)"
echo "Allocated CPUs: ${CPUS_PER_JOB}"
echo "Allocated Memory: ${MEM_PER_JOB}"
echo "Input FASTQ: ${FASTQ_FILE}"

# Load STAR module
module purge
module load STAR/2.7.11b # Or your preferred/available version

# Define fixed paths again inside the script for clarity
STAR_INDEX_JOB="${STAR_INDEX}"
OUTPUT_DIR_JOB="${OUTPUT_DIR}" # Output for BAM file

# Define output prefix path for STAR
OUT_PREFIX="\${OUTPUT_DIR_JOB}/${SAMPLE_ID}_"

# Run STAR alignment command
echo "Running STAR..."
STAR \\
    --genomeDir "\${STAR_INDEX_JOB}" \\
    --readFilesIn "${FASTQ_FILE}" \\
    --readFilesCommand zcat \\
    --runThreadN ${CPUS_PER_JOB} \\
    --outFileNamePrefix "\${OUT_PREFIX}" \\
    --outSAMtype BAM SortedByCoordinate \\
    --twopassMode Basic \\
    --limitBAMsortRAM ${SORT_RAM_LIMIT}

# Check STAR exit status
if [ \$? -eq 0 ]; then
    echo "STAR alignment completed successfully for sample ${SAMPLE_ID}"
else
    echo "STAR alignment failed for sample ${SAMPLE_ID}"
    # Optional: Add cleanup for failed job's partial output if needed
    # rm -f "\${OUT_PREFIX}"*
    exit 1
fi

echo "STAR processing complete for ${SAMPLE_ID}."

EOF

    # Make the generated script executable
    chmod +x "${SLURM_SCRIPT_NAME}"

done

echo "Finished generating scripts in ${SCRIPT_DIR}"
echo "To submit all jobs, you can run:"
echo "cd ${SCRIPT_DIR}"
echo "for script in *.sh; do sbatch \"\$script\"; done"

# /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_STAR_alignment_generate_star_scripts_for_part2_files.sh