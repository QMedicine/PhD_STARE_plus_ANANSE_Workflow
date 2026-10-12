#!/bin/bash
#SBATCH --job-name=run_manual_tfactivity
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/run_manual_tfactivity_logs/run_manual_tfactivity_res_100000_w_200000_m_2000000=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/run_manual_tfactivity_logs/run_manual_tfactivity_res_100000_w_200000_m_2000000=%x_%j.err

set -euo pipefail

module purge
module load miniconda3

# Define the conda function for this shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python3

# load R
module load R/4.5.1
export R_LIBS_USER="/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
CORES="${SLURM_CPUS_PER_TASK:-6}"

### ====== EDIT THESE ======
PROJ_DIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/Run_manual_tfactivity_DYNAMITE_ranking
RES=100000
RES_SHORT="100kb"
W_VAL=200000
M_VAL=2000000
WORK_DIR=${PROJ_DIR}/run_manual_tf_activity_results_res_${RES}_w_${W_VAL}_m_${M_VAL}

# Condition labels (used in filenames)
COND_A="Young_LTHSC"
COND_B="Old_LTHSC"
ASSAY="ATAC-seq"

COND_A_ALT="A_Young_HSC"
COND_B_ALT="B_Old_HSC"

TF_GENE_AFF_DIR_OLD=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Old_LTHSC/Gene_TF_matrices
TF_GENE_AFF_DIR_YOUNG=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Young_LTHSC/Gene_TF_matrices

# Glob (or space-separated list) of STARE TF–gene matrices per condition (your *_TF_Gene_Affinities_*.txt.gz)
COND_A_FILES=( "${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep1.txt.gz" \
"${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep2.txt.gz" )

COND_B_FILES=( "${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep1.txt.gz" \
"${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep2.txt.gz" )



set -x  # turn on debug printing (optional)
printf "YNG files:\n%s\n" "${COND_A_FILES[@]}"
printf "OLD files:\n%s\n" "${COND_B_FILES[@]}"
for f in "${COND_A_FILES[@]}" "${COND_B_FILES[@]}"; do ls -lh "$f"; done

# DESeq2 results (TSV with rownames=gene IDs, col 'log2FoldChange')
DESEQ2_TSV="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/${COND_A_ALT}:${COND_B_ALT}.deseq2.results.tsv"

# Optional: gene ID→symbol map (two columns: gene_id  gene_name).
# If you don't have one, provide a GTF and we'll build it:
ID_SYMBOL_MAP="https://ondemand.hpc.wehi.edu.au/pun/sys/dashboard/files/fs//vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_expanded_motif_dimers/00_prepare_genome/03_id_symbol_map/id_symbol_map.txt"     # or leave empty
GTF="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf"     # used only if ID_SYMBOL_MAP=""

### ========================

mkdir -p ${WORK_DIR}

# 0a) Build id_symbol_map.txt if not supplied
if [[ ! -s "${ID_SYMBOL_MAP:-}" ]]; then
  echo "[INFO] Building id_symbol_map.txt from ${GTF}"
  python3 ${PROJ_DIR}/00_make_id_symbol_map.py --gtf "$GTF" --out ${WORK_DIR}/id_symbol_map.txt
  ID_SYMBOL_MAP="${WORK_DIR}/id_symbol_map.txt"
fi

# 1) Combine replicates by mean (per condition)
python3 ${PROJ_DIR}/01_affinity_mean.py --method mean --inputs "${COND_A_FILES[@]}" --out "${WORK_DIR}/${COND_A}_${ASSAY}.tsv"
python3 ${PROJ_DIR}/01_affinity_mean.py --method mean --inputs "${COND_B_FILES[@]}" --out "${WORK_DIR}/${COND_B}_${ASSAY}.tsv"

# 2) Aggregate synonyms (genes & TFs)
python3 ${PROJ_DIR}/02_aggregate_synonyms.py \
  --in "${WORK_DIR}/${COND_A}_${ASSAY}.tsv" \
  --idmap "${ID_SYMBOL_MAP}" \
  --agg max \
  --out "${WORK_DIR}/${COND_A}_${ASSAY}.agg_affinities.tsv"

python3 ${PROJ_DIR}/02_aggregate_synonyms.py \
  --in "${WORK_DIR}/${COND_B}_${ASSAY}.tsv" \
  --idmap "${ID_SYMBOL_MAP}" \
  --agg max \
  --out "${WORK_DIR}/${COND_B}_${ASSAY}.agg_affinities.tsv"

# 3) Affinity ratio (B/A) with NA/Inf dropping & rounding (as in tfactivity)
python3 ${PROJ_DIR}/03_affinity_ratio.py \
  --a "${WORK_DIR}/${COND_A}_${ASSAY}.agg_affinities.tsv" \
  --b "${WORK_DIR}/${COND_B}_${ASSAY}.agg_affinities.tsv" \
  --out "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tsv"

# 4) DYNAMITE preprocess (merge with DE sign into 'Expression' 0/1)
python3 ${PROJ_DIR}/04_dynamite_preprocess.py \
  --ratio "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tsv" \
  --deseq "${DESEQ2_TSV}" \
  --out "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.preprocessed.tsv"

# 5) Run DYNAMITE (elastic-net logistic CV) and move coefficients
bash ${PROJ_DIR}/05_run_dynamite.sh \
  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.preprocessed.tsv" \
  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}_dynamite_regression_coefficients.txt" \
  3 6 0.1 "${CORES}" false

# 6) Filter coefficients by |value| ≥ 0.1 (same awk logic)
python3 ${PROJ_DIR}/06_filter_coefficients.py \
  --in  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}_dynamite_regression_coefficients.txt" \
  --thr 0.1 \
  --out "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.filtered.tsv"

# 7) TF–TG score = |log2FC| × affinity × |coef|
python3 ${PROJ_DIR}/07_tf_tg_score.py \
  --deseq  "${DESEQ2_TSV}" \
  --aff    "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tsv" \
  --coefs  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.filtered.tsv" \
  --out    "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.score.tsv"

# 8) Rank TFs (MWU vs background + DCG) and gene-wise DCG maps
python3 ${PROJ_DIR}/08_create_ranking.py \
  --scores "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.score.tsv" \
  --tfout  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.tsv" \
  --tgout  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tg_ranking.tsv"

echo "[DONE] Key outputs:"
echo " - TF ranking: ${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.tsv"
echo " - TG (per-TF gene DCG): ${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tg_ranking.tsv"
echo " - DYNAMITE coefficients (filtered): ${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.filtered.tsv"

# 9)
python3 ${PROJ_DIR}/09_build_tf_mapping_and_annotate.py \
  --inputs \
    "${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep1.txt.gz" \
    "${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep2.txt.gz" \
    "${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep1.txt.gz" \
    "${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep2.txt.gz" \
  --motifs "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_expanded_motif_dimers/02_motifs/05_transfac/motifs.converted.transfac" \
  --expanded_map_csv "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/jaspar_expanded_map.csv" \
  --tf_ranking "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.tsv" \
  --out_map "${WORK_DIR}/tf_mapping.tsv" \
  --out_tf_rank_annot "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.annot.tsv"

# 10) Calculate ABC score ratios based on STARE output
ABC_OLD_DIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Old_LTHSC/ABC_output
ABC_YOUNG_DIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Young_LTHSC/ABC_output

ABC_OLD_FILES=(
  "${ABC_OLD_DIR}/Old_LTHSC_ABCpp_scoredInteractions_Old_rep1.txt.gz"
  "${ABC_OLD_DIR}/Old_LTHSC_ABCpp_scoredInteractions_Old_rep2.txt.gz"
)
ABC_YOUNG_FILES=(
  "${ABC_YOUNG_DIR}/Young_LTHSC_ABCpp_scoredInteractions_Young_rep1.txt.gz"
  "${ABC_YOUNG_DIR}/Young_LTHSC_ABCpp_scoredInteractions_Young_rep2.txt.gz"
)

for f in "${ABC_OLD_FILES[@]}" "${ABC_YOUNG_FILES[@]}"; do
  [[ -s "$f" ]] || { echo "[ERR] Missing: $f" >&2; exit 2; }
done

python3 ${PROJ_DIR}/10_compute_gene_gabc_ratio.py \
  --abc_old   "${ABC_OLD_FILES[@]}" \
  --abc_young "${ABC_YOUNG_FILES[@]}" \
  --agg mean \
  --gene-col "Gene Name" \
  --score-col "ABC-Score" \
  --out "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.log2_gabc_ratio.tsv"
  
# 11) Summarize all results into a table
python3 ${PROJ_DIR}/11_summarize_top_tfs.py \
  --tf_rank "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.annot.tsv" \
  --score   "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.score.tsv" \
  --coef    "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.filtered.tsv" \
  --ratio   "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tsv" \
  --deseq   "${DESEQ2_TSV}" \
  --gabc    "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.log2_gabc_ratio.tsv" \
  --topN 20 --topM 25 --delta 0.25 \
  --tf_tpm /vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window/01_counts/04_filtered_tfs/TFs.tpm_filtered.tsv \
  --idmap   "${ID_SYMBOL_MAP}" \
  --expanded-map "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/motifs/jaspar_expanded_map.csv" \
  --tf_map "${WORK_DIR}/tf_mapping.tsv" \
  --out     "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.topTF_topGenes.summary.with_gabc.tsv"

echo "finished summarizing all results into a table"

# 12)
python3 ${PROJ_DIR}/12_plot_global_tf_ranking.py \
  --tf_rank_annot "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.annot.tsv" \
  --coefs         "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.filtered.tsv" \
  --summary "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.topTF_topGenes.summary.with_gabc.tsv" \
  --topN 20 \
  --out "${WORK_DIR}/fig_global_tf_ranking_top20.pdf"

# 13) 
python3 ${PROJ_DIR}/13_plot_tf_gene_heatmap.py \
  --summary        "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.topTF_topGenes.summary.with_gabc.tsv" \
  --tf_rank_annot  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.annot.tsv" \
  --topT 10 --topG 25 \
  --out "${WORK_DIR}/fig_tf_gene_heatmap_score_T10xG20.pdf"

# 14)
mkdir -p "${WORK_DIR}/tf_panels"
python3 ${PROJ_DIR}/14_plot_per_tf_panels.py \
  --summary        "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.topTF_topGenes.summary.with_gabc.tsv" \
  --tf_rank_annot  "${WORK_DIR}/${COND_A}:${COND_B}_${ASSAY}.tf_ranking.annot.tsv" \
  --tf_tpm         "/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window/01_counts/04_filtered_tfs/TFs.tpm_filtered.tsv" \
  --deseq          "${DESEQ2_TSV}" \
  --idmap          "${ID_SYMBOL_MAP}" \
  --topTF 6 --topG 25 \
  --label-top 10 --label-metric score \
  --outdir "${WORK_DIR}/tf_panels"

echo "plotting finished"