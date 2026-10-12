#!/bin/bash
#SBATCH --job-name=index_bam_files
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/index_bam_log/index_bam_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/index_bam_log/index_bam_%A_%a.err
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=1:00:00

# Load samtools module
module purge
module load samtools

# Check if samtools is available in the environment
echo "Checking samtools availability..."
which samtools

# Define directories
BAM_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams"

# Get the list of BAM files
BAM_FILES=(${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

# Select the current BAM file based on the job array index
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the sample ID
SAMPLE_ID=$(basename ${BAM_FILE} "_Aligned.sortedByCoord.out.bam")

# Debug information
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Indexing BAM file: ${BAM_FILE}"
echo "Sample ID: ${SAMPLE_ID}"

# Index the BAM file
samtools index ${BAM_FILE}

# Check if samtools index completed successfully
if [ $? -eq 0 ]; then
    echo "Indexing completed successfully for sample ${SAMPLE_ID}"
else
    echo "Indexing failed for sample ${SAMPLE_ID}"
    exit 1
fi

echo "Indexing of sample ${SAMPLE_ID} is complete."


# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_index_bam.sh

