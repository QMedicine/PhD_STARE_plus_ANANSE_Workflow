#!/bin/bash
module load nextflow/25.04.2

# ---------- paths ----------
PROJ_DIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity
SCRATCH=/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_nfcore_tfactivity
OUTDIR=${PROJ_DIR}/tfactivity_results/HSC_50kbp_window_JASPAR2024_all_vert_motifs_Young_over_Old
mkdir -p "$OUTDIR"
WORKDIR=${SCRATCH}/work_50kbp_window_JASPAR2024_all_vert_motifs_Young_over_Old

# ---------- container caches ----------
export NXF_SINGULARITY_CACHEDIR=/vast/projects/Immune_Aging/apptainer_singularity_cache
export APPTAINER_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
export SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
export APPTAINER_TMPDIR=/vast/scratch/users/stacpoole.q/apptainer_singularity_cache
mkdir -p "$APPTAINER_TMPDIR/tmp"

# ---------- JVM & bashrc ----------
source ~/.bashrc
export NXF_OPTS='-Xms1g -Xmx4g'

# ---------- run ----------
nextflow -log $PROJ_DIR/nf_logs/nf_tfactivity_50kbp_window_JASPAR2024_all_vert_motifs_Young_over_Old.log run nf-core/tfactivity -r dev \
       -profile wehi \
       -c  $PROJ_DIR/tfactivity_timeout.config \
       -resume \
       -w "${WORKDIR}" \
       --input        "${PROJ_DIR}/samplesheet_HSC_Young_over_Old.csv" \
       --fasta        /vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa \
       --gtf          /vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf \
       --blacklist    /vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed \
       --motifs     "${PROJ_DIR}/motifs/JASPAR2024_all_vertebrates_motifs.jaspar" \
       --counts       "/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv" \
       --counts_design "${PROJ_DIR}/design_matrix_HSC_Young_over_Old.csv" \
       --min_tpm 0.1 \
       --min_tpm_tf 0.1 \
       --duplicate_motifs merge \
       --outdir       "${OUTDIR}" \
       --window_size 50000 \
       --max_cpus 32 --max_memory 380.GB --max_time 48.h --skip_sneep \
       -bg --email stacpoole.q@wehi.edu.au