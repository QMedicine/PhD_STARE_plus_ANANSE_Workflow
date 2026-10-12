#!/bin/bash
#SBATCH --job-name=wig_to_bigwig_conversion
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/ucsctools_logs/wig_to_bigwig_conversion_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/ucsctools_logs/wig_to_bigwig_conversion_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00

# Load the UCSC tools module or set the path to wigToBigWig
module load ucsc-tools

# Define the directory containing the WIG files (with .bw extension)
WIG_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/wiggletools_combined_coverage_tracks"

# Path to the chromosome sizes file
CHROM_SIZES="${WIG_DIR}/chrom.sizes.txt"

# Check if the chromosome sizes file exists, generate it if not
if [ ! -f "${CHROM_SIZES}" ]; then
    echo "Chromosome sizes file not found at ${CHROM_SIZES}"
    echo "Generating chromosome sizes file..."
    # Define the path to the genome fasta file
    GENOME_FASTA="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa"
    if [ ! -f "${GENOME_FASTA}" ]; then
        echo "Genome FASTA file not found at ${GENOME_FASTA}"
        exit 1
    fi
    samtools faidx "${GENOME_FASTA}"
    cut -f1,2 "${GENOME_FASTA}.fai" > "${CHROM_SIZES}"
    if [ ! -f "${CHROM_SIZES}" ]; then
        echo "Failed to create chromosome sizes file."
        exit 1
    fi
fi

# Check if wigToBigWig is available
if ! command -v wigToBigWig &> /dev/null; then
    echo "wigToBigWig could not be found. Please ensure it is installed and in your PATH."
    exit 1
fi

# Loop over all .bw files in the directory (which are actually WIG files)
for WIG_FILE in "${WIG_DIR}"/*.bw; do
    # Check if there are any files matching the pattern
    if [ ! -e "${WIG_FILE}" ]; then
        echo "No .bw files found in ${WIG_DIR}"
        exit 1
    fi
    
    # Get the base name of the file
    BASE_NAME=$(basename "${WIG_FILE}")
    # Remove the .bw extension
    FILE_NAME="${BASE_NAME%.bw}"
    # Define the output file name
    OUTPUT_FILE="${WIG_DIR}/${FILE_NAME}_converted.bw"
    
    echo "Converting ${WIG_FILE} to ${OUTPUT_FILE}"
    
    # Run wigToBigWig
    wigToBigWig "${WIG_FILE}" "${CHROM_SIZES}" "${OUTPUT_FILE}"
    
    # Check if the output file was created
    if [ -f "${OUTPUT_FILE}" ]; then
        echo "Successfully created ${OUTPUT_FILE}"
    else
        echo "Failed to create ${OUTPUT_FILE}"
        exit 1
    fi
done

echo "Conversion completed successfully."

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_convert_wig_to_bigwig.sh