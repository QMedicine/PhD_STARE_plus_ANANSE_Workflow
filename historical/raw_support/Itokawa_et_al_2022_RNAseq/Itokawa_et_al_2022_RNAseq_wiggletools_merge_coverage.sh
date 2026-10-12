#!/bin/bash
#SBATCH --job-name=wiggletools_average_coverage_tracks
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/wiggletoole_logs/average_coverage_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/wiggletoole_logs/average_coverage_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00

# Define environment path (assuming the environment already exists)
ENV_PATH="/home/users/allstaff/stacpoole.q/conda_envs/wiggletools_env"

# Define the path to wiggletools
WIGGLETOOLS_BIN="${ENV_PATH}/bin/wiggletools"

# Check if wiggletools exists
if [ ! -x "${WIGGLETOOLS_BIN}" ]; then
    echo "wiggletools not found at ${WIGGLETOOLS_BIN}"
    exit 1
fi

# Load samtools module for generating chrom.sizes
module load samtools

# Define directories
COVERAGE_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/deepTools_coverage_tracks"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/wiggletools_combined_coverage_tracks"
GENOME_FASTA="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa"
CHROM_SIZES="${OUTPUT_DIR}/chrom.sizes.txt"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Generate chrom.sizes file if it doesn't exist
if [ ! -f ${CHROM_SIZES} ]; then
    echo "Generating chrom.sizes file..."
    samtools faidx ${GENOME_FASTA}
    cut -f1,2 ${GENOME_FASTA}.fai > ${CHROM_SIZES}
fi

# Create arrays of BigWig files for each condition
YOUNG_BWS=(${COVERAGE_DIR}/*Young_HSC*_deepTools.bw)
OLD_BWS=(${COVERAGE_DIR}/*Old_HSC*_deepTools.bw)

# Check if BigWig files exist for each condition
if [ ${#YOUNG_BWS[@]} -eq 0 ]; then
    echo "No 'Young_HSC' replicate BigWig files found in ${COVERAGE_DIR}"
    exit 1
fi

if [ ${#OLD_BWS[@]} -eq 0 ]; then
    echo "No 'Old_HSC' replicate BigWig files found in ${COVERAGE_DIR}"
    exit 1
fi

# Debugging: List the BigWig files
echo "Young_HSC BigWig files:"
printf '%s\n' "${YOUNG_BWS[@]}"

echo "Old_HSC BigWig files:"
printf '%s\n' "${OLD_BWS[@]}"

# Define output files
YOUNG_OUTPUT="${OUTPUT_DIR}/Young_HSC_combined.bw"
OLD_OUTPUT="${OUTPUT_DIR}/Old_HSC_combined.bw"

# Remove existing output files if they exist
if [ -f "${YOUNG_OUTPUT}" ]; then
    echo "Deleting existing file ${YOUNG_OUTPUT}"
    rm "${YOUNG_OUTPUT}"
fi

if [ -f "${OLD_OUTPUT}" ]; then
    echo "Deleting existing file ${OLD_OUTPUT}"
    rm "${OLD_OUTPUT}"
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_HSC' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_OUTPUT}" mean "${YOUNG_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_HSC' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_OUTPUT}" mean "${OLD_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_OUTPUT}" ]; then
    echo "Failed to create ${OLD_OUTPUT}"
    exit 1
fi

echo "Averaging of coverage tracks completed successfully."


# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_wiggletools_merge_coverage.sh