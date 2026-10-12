#!/bin/bash
#SBATCH --job-name=STAR_alignment
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_alignment_log/STAR_alignment_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_alignment_log/STAR_alignment_%A_%a.err
#SBATCH --array=0-3  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time=24:00:00

# Load STAR module
module purge
module load STAR

# Check if STAR is available in the environment
echo "Checking STAR availability..."
which STAR

# Print the STAR version
STAR --version

# Define directories
TRIMMED_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed"
STAR_INDEX="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/STAR_index_60bp"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams"
METRICS_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_alignment_metrics"

# Create output and metrics directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${METRICS_DIR}

# Get the list of FASTQ files
FASTQ_FILES=(${TRIMMED_DIR}/*_trimmed.fastq.gz)

# Select the current FASTQ file based on the job array index
FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the sample ID
SAMPLE_ID=$(basename ${FASTQ_FILE} "_trimmed.fastq.gz")

# Debug information
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing file: ${FASTQ_FILE}"
echo "Sample ID: ${SAMPLE_ID}"

# Run STAR alignment
STAR \
    --genomeDir ${STAR_INDEX} \
    --readFilesIn ${FASTQ_FILE} \
    --readFilesCommand zcat \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE_ID}_ \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --limitBAMsortRAM 20000000000

# Check if STAR completed successfully
if [ $? -eq 0 ]; then
    echo "STAR alignment completed successfully for sample ${SAMPLE_ID}"
else
    echo "STAR alignment failed for sample ${SAMPLE_ID}"
    exit 1
fi

echo "Processing of sample ${SAMPLE_ID} is complete."



# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_STAR_alignment.sh


