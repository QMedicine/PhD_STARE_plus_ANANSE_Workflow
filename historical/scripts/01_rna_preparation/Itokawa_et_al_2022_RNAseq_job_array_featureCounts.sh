#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/featureCounts_logs/featureCounts_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/featureCounts_logs/featureCounts_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Load the required module
module purge
module load subread

# Define variables
ANNOTATION="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf"
BAM_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams"
OUTPUT_DIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/featureCounts"
mkdir -p ${OUTPUT_DIR}

# Generate a list of BAM files
BAM_FILES=(${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

# Run featureCounts
featureCounts -T 8 -a ${ANNOTATION} --extraAttributes gene_name,gene_biotype -o ${OUTPUT_DIR}/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt ${BAM_FILES[@]}

# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_featureCounts.sh