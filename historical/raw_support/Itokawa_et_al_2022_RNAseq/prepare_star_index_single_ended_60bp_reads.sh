#!/bin/bash

# Define paths and parameters
GENOME_DIR="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113"

# Path to the genome FASTA file (Ensembl GRCm39 primary assembly)
FASTA_FILE="${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# Path to the gene annotation GTF file (Ensembl release 113)
GTF_FILE="${GENOME_DIR}/Mus_musculus.GRCm39.113.gtf"

# Directory where the STAR genome index will be generated
STAR_INDEX_DIR="${GENOME_DIR}/STAR_index_60bp"

# Number of threads to use during genome indexing
THREADS=8

# Read length of the RNA-seq data
READ_LENGTH=60  # Set this based on the RNA-seq read length

# Calculate sjdbOverhang parameter (ideal value is READ_LENGTH - 1)
SJDB_OVERHANG=$((READ_LENGTH - 1))

# Create the STAR index directory if it doesn't exist
mkdir -p ${STAR_INDEX_DIR}

module load STAR

# Run STAR to generate the genome index
STAR --runMode genomeGenerate \
     --genomeDir ${STAR_INDEX_DIR} \
     --genomeFastaFiles ${FASTA_FILE} \
     --sjdbGTFfile ${GTF_FILE} \
     --sjdbOverhang ${SJDB_OVERHANG} \
     --runThreadN ${THREADS}

# Check if STAR completed successfully and print a message
if [ $? -eq 0 ]; then
    echo "STAR genome index successfully created in: ${STAR_INDEX_DIR}"
else
    echo "STAR genome index generation failed."
fi

#chmod +x /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/prepare_star_index_single_ended_60bp_reads.sh