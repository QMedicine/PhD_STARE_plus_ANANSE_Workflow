#!/bin/bash

# Set paths and parameters
TRANSCRIPTOME_DIR="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113"
TRANSCRIPTOME_FASTA="${TRANSCRIPTOME_DIR}/Mus_musculus.GRCm39.cdna.all.fa"
SALMON_INDEX_DIR="${TRANSCRIPTOME_DIR}/Salmon_index"
THREADS=8  # Adjust based on your available resources

# Create the Salmon index directory if it doesn't exist
mkdir -p ${SALMON_INDEX_DIR}

module purge
module load salmon

# Build the Salmon index
salmon index -t ${TRANSCRIPTOME_FASTA} \
             -i ${SALMON_INDEX_DIR} \
             -k 31 \
             --threads ${THREADS}

# Print a message upon completion
if [ $? -eq 0 ]; then
    echo "Salmon index successfully created in: ${SALMON_INDEX_DIR}"
else
    echo "Salmon index generation failed."
fi

# chmod +x /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/prepare_salmon_index_single_ended_60bp_reads.sh