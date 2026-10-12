#!/bin/bash
#SBATCH --job-name=deepTools_coverage_tracks
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/deepTools_coverage_log/deepTools_generate_coverage_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/deepTools_coverage_log/deepTools_generate_coverage_%A_%a.err
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=8:00:00

# Load deepTools module
module purge
module load deeptools

# Define directories
BAM_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams"
COVERAGE_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/deepTools_coverage_tracks"

# Create coverage directory if it doesn't exist
mkdir -p ${COVERAGE_DIR}

# Get the list of BAM files
BAM_FILES=(${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

# Select the current BAM file based on the job array index
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the sample ID
SAMPLE_ID=$(basename ${BAM_FILE} "_Aligned.sortedByCoord.out.bam")

# Debug information
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing BAM file: ${BAM_FILE}"
echo "Sample ID: ${SAMPLE_ID}"

# Run bamCoverage
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${COVERAGE_DIR}/${SAMPLE_ID}_deepTools.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 4

# Check if bamCoverage completed successfully
if [ $? -eq 0 ]; then
    echo "Coverage track generated successfully for sample ${SAMPLE_ID}"
else
    echo "Coverage track generation failed for sample ${SAMPLE_ID}"
    exit 1
fi

echo "Processing of sample ${SAMPLE_ID} is complete."

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_deepTools_coverage.sh



