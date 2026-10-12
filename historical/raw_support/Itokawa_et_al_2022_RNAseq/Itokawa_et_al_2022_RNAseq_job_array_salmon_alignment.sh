#!/bin/bash
#SBATCH --job-name=STAR_alignment
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/salmon_alignment_log/salmon_alignment_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/salmon_alignment_log/salmon_alignment_%A_%a.err
#SBATCH --array=0-3  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# Load Salmon module
module purge
module load salmon

# Check if Salmon is available in the environment
echo "Checking Salmon availability..."
which salmon

# Print the Salmon version
salmon --version

# Define directories
TRIMMED_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/trimmed"
SALMON_INDEX="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Salmon_index"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/salmon_quant"
METRICS_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/salmon_metrics"

# Create output and metrics directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${METRICS_DIR}

# Get the list of FASTQ files
FASTQ_FILES=(${TRIMMED_DIR}/*_trimmed.fastq.gz)

# Select the current FASTQ file based on the job array index
FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the sample ID
SAMPLE_ID=$(basename ${FASTQ_FILE} "_trimmed.fastq.gz")

# Define output directories for the sample
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
SAMPLE_METRICS_FILE="${METRICS_DIR}/${SAMPLE_ID}_metrics.txt"

# Debug information
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing file: ${FASTQ_FILE}"
echo "Sample ID: ${SAMPLE_ID}"
echo "Salmon output directory: ${SAMPLE_OUTPUT_DIR}"

# Run Salmon quantification
salmon quant \
    -i ${SALMON_INDEX} \
    --libType SR \
    -r ${FASTQ_FILE} \
    --seqBias \
    --gcBias \
    --posBias \
    --writeUnmappedNames \
    --threads 8 \
    -o ${SAMPLE_OUTPUT_DIR}

# Capture metrics
if [ $? -eq 0 ]; then
    grep -E "Mapping rate|Fragmentation|Duplicate" ${SAMPLE_OUTPUT_DIR}/logs/salmon_quant.log > ${SAMPLE_METRICS_FILE}
    echo "Salmon quantification completed successfully for ${SAMPLE_ID}."
    echo "Metrics saved to: ${SAMPLE_METRICS_FILE}"
else
    echo "Salmon quantification failed for ${SAMPLE_ID}."
fi

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_salmon_alignment.sh