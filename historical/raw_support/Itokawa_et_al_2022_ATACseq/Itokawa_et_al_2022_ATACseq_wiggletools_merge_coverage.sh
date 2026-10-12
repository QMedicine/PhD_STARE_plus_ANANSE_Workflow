#!/bin/bash
#SBATCH --job-name=wiggletools_average_coverage_tracks
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/wiggletools_logs/average_coverage_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/wiggletools_logs/average_coverage_%j.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

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
COVERAGE_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/deepTools_coverage_tracks"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/wiggletools_combined_coverage_tracks"
GENOME_FASTA="/vast/projects/Immune_Aging/Reference_genomes/mm39/Bowtie2/Mus_musculus.GRCm39.dna.primary_assembly.fa"
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
YOUNG_HSC_BWS=(${COVERAGE_DIR}/*Young_HSC*_deepTools.bw)
OLD_HSC_BWS=(${COVERAGE_DIR}/*Old_HSC*_deepTools.bw)

YOUNG_MPP1_BWS=(${COVERAGE_DIR}/*Young_MPP1*_deepTools.bw)
OLD_MPP1_BWS=(${COVERAGE_DIR}/*Old_MPP1*_deepTools.bw)

YOUNG_MPP2_BWS=(${COVERAGE_DIR}/*Young_MPP2*_deepTools.bw)
OLD_MPP2_BWS=(${COVERAGE_DIR}/*Old_MPP2*_deepTools.bw)

YOUNG_MPP3_BWS=(${COVERAGE_DIR}/*Young_MPP3*_deepTools.bw)
OLD_MPP3_BWS=(${COVERAGE_DIR}/*Old_MPP3*_deepTools.bw)

YOUNG_MPP4_BWS=(${COVERAGE_DIR}/*Young_MPP4*_deepTools.bw)
OLD_MPP4_BWS=(${COVERAGE_DIR}/*Old_MPP4*_deepTools.bw)

YOUNG_GMP_BWS=(${COVERAGE_DIR}/*Young_GMP*_deepTools.bw)
OLD_GMP_BWS=(${COVERAGE_DIR}/*Old_GMP*_deepTools.bw)

YOUNG_MEP_BWS=(${COVERAGE_DIR}/*Young_MEP*_deepTools.bw)
OLD_MEP_BWS=(${COVERAGE_DIR}/*Old_MEP*_deepTools.bw)

YOUNG_CLP_BWS=(${COVERAGE_DIR}/*Young_CLP*_deepTools.bw)
OLD_CLP_BWS=(${COVERAGE_DIR}/*Old_CLP*_deepTools.bw)

# Check if BigWig files exist for each condition
if [ ${#YOUNG_HSC_BWS[@]} -eq 0 ]; then
    echo "No 'Young_HSC' replicate BigWig files found in ${COVERAGE_DIR}"
    exit 1
fi

if [ ${#OLD_HSC_BWS[@]} -eq 0 ]; then
    echo "No 'Old_HSC' replicate BigWig files found in ${COVERAGE_DIR}"
    exit 1
fi

# Debugging: List the BigWig files
echo "Young_HSC BigWig files:"
printf '%s\n' "${YOUNG_HSC_BWS[@]}"

echo "Old_HSC BigWig files:"
printf '%s\n' "${OLD_HSC_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_MPP1 BigWig files:"
printf '%s\n' "${YOUNG_MPP1_BWS[@]}"

echo "Old_MPP1 BigWig files:"
printf '%s\n' "${OLD_MPP1_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_MPP2 BigWig files:"
printf '%s\n' "${YOUNG_MPP2_BWS[@]}"

echo "Old_MPP2 BigWig files:"
printf '%s\n' "${OLD_MPP2_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_MPP3 BigWig files:"
printf '%s\n' "${YOUNG_MPP3_BWS[@]}"

echo "Old_MPP3 BigWig files:"
printf '%s\n' "${OLD_MPP3_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_MPP4 BigWig files:"
printf '%s\n' "${YOUNG_MPP4_BWS[@]}"

echo "Old_MPP4 BigWig files:"
printf '%s\n' "${OLD_MPP4_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_GMP BigWig files:"
printf '%s\n' "${YOUNG_GMP_BWS[@]}"

echo "Old_GMP BigWig files:"
printf '%s\n' "${OLD_GMP_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_MEP BigWig files:"
printf '%s\n' "${YOUNG_MEP_BWS[@]}"

echo "Old_MEP BigWig files:"
printf '%s\n' "${OLD_MEP_BWS[@]}"

# Debugging: List the BigWig files
echo "Young_CLP BigWig files:"
printf '%s\n' "${YOUNG_CLP_BWS[@]}"

echo "Old_CLP BigWig files:"
printf '%s\n' "${OLD_CLP_BWS[@]}"

# Define output files
YOUNG_HSC_OUTPUT="${OUTPUT_DIR}/Young_HSC_combined.bw"
OLD_HSC_OUTPUT="${OUTPUT_DIR}/Old_HSC_combined.bw"

YOUNG_MPP1_OUTPUT="${OUTPUT_DIR}/Young_MPP1_combined.bw"
OLD_MPP1_OUTPUT="${OUTPUT_DIR}/Old_MPP1_combined.bw"

YOUNG_MPP2_OUTPUT="${OUTPUT_DIR}/Young_MPP2_combined.bw"
OLD_MPP2_OUTPUT="${OUTPUT_DIR}/Old_MPP2_combined.bw"

YOUNG_MPP3_OUTPUT="${OUTPUT_DIR}/Young_MPP3_combined.bw"
OLD_MPP3_OUTPUT="${OUTPUT_DIR}/Old_MPP3_combined.bw"

YOUNG_MPP4_OUTPUT="${OUTPUT_DIR}/Young_MPP4_combined.bw"
OLD_MPP4_OUTPUT="${OUTPUT_DIR}/Old_MPP4_combined.bw"

YOUNG_GMP_OUTPUT="${OUTPUT_DIR}/Young_GMP_combined.bw"
OLD_GMP_OUTPUT="${OUTPUT_DIR}/Old_GMP_combined.bw"

YOUNG_MEP_OUTPUT="${OUTPUT_DIR}/Young_MEP_combined.bw"
OLD_MEP_OUTPUT="${OUTPUT_DIR}/Old_MEP_combined.bw"

YOUNG_CLP_OUTPUT="${OUTPUT_DIR}/Young_CLP_combined.bw"
OLD_CLP_OUTPUT="${OUTPUT_DIR}/Old_CLP_combined.bw"

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_HSC' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_HSC_OUTPUT}" mean "${YOUNG_HSC_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_HSC_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_HSC_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_HSC' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_HSC_OUTPUT}" mean "${OLD_HSC_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_HSC_OUTPUT}" ]; then
    echo "Failed to create ${OLD_HSC_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_MPP1' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_MPP1_OUTPUT}" mean "${YOUNG_MPP1_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_MPP1_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_MPP1_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_MPP1' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_MPP1_OUTPUT}" mean "${OLD_MPP1_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_MPP1_OUTPUT}" ]; then
    echo "Failed to create ${OLD_MPP1_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_MPP2' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_MPP2_OUTPUT}" mean "${YOUNG_MPP2_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_MPP2_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_MPP2_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_MPP2' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_MPP2_OUTPUT}" mean "${OLD_MPP2_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_MPP2_OUTPUT}" ]; then
    echo "Failed to create ${OLD_MPP2_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_MPP3' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_MPP3_OUTPUT}" mean "${YOUNG_MPP3_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_MPP3_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_MPP3_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_MPP3' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_MPP3_OUTPUT}" mean "${OLD_MPP3_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_MPP3_OUTPUT}" ]; then
    echo "Failed to create ${OLD_MPP3_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_MPP4' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_MPP4_OUTPUT}" mean "${YOUNG_MPP4_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_MPP4_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_MPP4_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_MPP4' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_MPP4_OUTPUT}" mean "${OLD_MPP4_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_MPP4_OUTPUT}" ]; then
    echo "Failed to create ${OLD_MPP4_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_GMP' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_GMP_OUTPUT}" mean "${YOUNG_GMP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_GMP_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_GMP_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_GMP' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_GMP_OUTPUT}" mean "${OLD_GMP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_GMP_OUTPUT}" ]; then
    echo "Failed to create ${OLD_GMP_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_MEP' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_MEP_OUTPUT}" mean "${YOUNG_MEP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_MEP_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_MEP_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_MEP' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_MEP_OUTPUT}" mean "${OLD_MEP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_MEP_OUTPUT}" ]; then
    echo "Failed to create ${OLD_MEP_OUTPUT}"
    exit 1
fi

# Average coverage for the young condition and write directly to BigWig
echo "Averaging coverage for the 'Young_CLP' condition..."
"${WIGGLETOOLS_BIN}" write "${YOUNG_CLP_OUTPUT}" mean "${YOUNG_CLP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${YOUNG_CLP_OUTPUT}" ]; then
    echo "Failed to create ${YOUNG_CLP_OUTPUT}"
    exit 1
fi

# Average coverage for the old condition and write directly to BigWig
echo "Averaging coverage for the 'Old_CLP' condition..."
"${WIGGLETOOLS_BIN}" write "${OLD_CLP_OUTPUT}" mean "${OLD_CLP_BWS[@]}"

# Check if the output file was created
if [ ! -f "${OLD_CLP_OUTPUT}" ]; then
    echo "Failed to create ${OLD_CLP_OUTPUT}"
    exit 1
fi

echo "Averaging of coverage tracks completed successfully."

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/Itokawa_et_al_2022_ATACseq_wiggletools_merge_coverage.sh