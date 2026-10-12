#!/bin/bash
#SBATCH --job-name=COsort_index_bam_files
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/COsort_index_bam_log/COsort_index_bam_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/COsort_index_bam_log/COsort_index_bam_%A_%a.err
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4            # Increased CPUs for multithreading
#SBATCH --mem=60G                    # Increased memory
#SBATCH --time=12:00:00

set -e  # Exit immediately if a command exits with a non-zero status

# Load samtools module
module purge
module load samtools

# Check if samtools is available in the environment
echo "Checking samtools availability..."
which samtools
samtools --version

# Define directories
BAM_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/filtered_bam_files"
SCRATCH_CO_SORT_BAM_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files"

# Create output directories if they don't exist
mkdir -p "${SCRATCH_CO_SORT_BAM_DIR}"

# Get the list of BAM files
BAM_FILES=(${BAM_DIR}/*_aligned_FMI_MD_sort.bam)
TOTAL_FILES=${#BAM_FILES[@]}
echo "Total number of BAM files found: ${TOTAL_FILES}"

# Check if SLURM_ARRAY_TASK_ID is within range
if [ ${SLURM_ARRAY_TASK_ID} -ge ${TOTAL_FILES} ]; then
    echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} exceeds number of BAM files ${TOTAL_FILES}. Exiting."
    exit 1
fi

# Select the current BAM file based on the job array index
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract the sample ID
SAMPLE_ID=$(basename "${BAM_FILE}" "_aligned_FMI_MD_sort.bam")

# Debug information
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sorting BAM file: ${BAM_FILE}"
echo "Sample ID: ${SAMPLE_ID}"

# Define the output sorted BAM file
SORTED_BAM_FILE="${SCRATCH_CO_SORT_BAM_DIR}/${SAMPLE_ID}_COsorted.bam"

# Check if the sorted BAM file already exists
if [ -e "${SORTED_BAM_FILE}" ]; then
    echo "Sorted BAM file ${SORTED_BAM_FILE} already exists. Skipping."
    exit 0
fi

# Sort the BAM file by coordinate
echo "Starting sorting for sample ${SAMPLE_ID}..."
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${SORTED_BAM_FILE}" "${BAM_FILE}"
echo "Sorting completed successfully for sample ${SAMPLE_ID}"

# Index the sorted BAM file
echo "Starting indexing for sample ${SAMPLE_ID}..."
samtools index -@ ${SLURM_CPUS_PER_TASK} "${SORTED_BAM_FILE}"
echo "Indexing completed successfully for sample ${SAMPLE_ID}"

echo "Processing of sample ${SAMPLE_ID} is complete."

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_Hi_Lo_HSC_COsort_index_bam.sh


