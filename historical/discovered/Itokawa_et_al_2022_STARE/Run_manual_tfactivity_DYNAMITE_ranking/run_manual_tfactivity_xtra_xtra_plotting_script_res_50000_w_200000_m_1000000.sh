#!/bin/bash
#SBATCH --job-name=run_manual_tfactivity_plot
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/run_manual_tfactivity_xtra_xtra_plotting_logs/run_manual_tfactivity_plots_res_50000_w_200000_m_1000000=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/run_manual_tfactivity_xtra_xtra_plotting_logs/run_manual_tfactivity_plots_res_50000_w_200000_m_1000000=%x_%j.err

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
RES=50000
RES_SHORT="50kb"
W_VAL=200000
M_VAL=1000000
WORK_DIR=${PROJ_DIR}/run_manual_tf_activity_results_res_${RES}_w_${W_VAL}_m_${M_VAL}

# Condition labels (used in filenames)
COND_A="Old_LTHSC"
COND_B="Young_LTHSC"
ASSAY="ATAC-seq"

COND_A_ALT="Old_HSC"
COND_B_ALT="Young_HSC"

TF_GENE_AFF_DIR_OLD=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Old_LTHSC/Gene_TF_matrices
TF_GENE_AFF_DIR_YOUNG=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES_SHORT}_w${W_VAL}_m${M_VAL}/Young_LTHSC/Gene_TF_matrices

# Glob (or space-separated list) of STARE TF–gene matrices per condition (your *_TF_Gene_Affinities_*.txt.gz)
COND_A_FILES=( "${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep1.txt.gz" \
"${TF_GENE_AFF_DIR_OLD}/Old_LTHSC_TF_Gene_Affinities_Old_rep2.txt.gz" )

COND_B_FILES=( "${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep1.txt.gz" \
"${TF_GENE_AFF_DIR_YOUNG}/Young_LTHSC_TF_Gene_Affinities_Young_rep2.txt.gz" )

set -x  # turn on debug printing (optional)
printf "OLD files:\n%s\n" "${COND_A_FILES[@]}"
printf "YNG files:\n%s\n" "${COND_B_FILES[@]}"
for f in "${COND_A_FILES[@]}" "${COND_B_FILES[@]}"; do ls -lh "$f"; done

# DESeq2 results (TSV with rownames=gene IDs, col 'log2FoldChange')
DESEQ2_TSV="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_1Mbp_window/01_counts/05_deseq2/Old_HSC:Young_HSC/${COND_A_ALT}:${COND_B_ALT}.deseq2.results.tsv"

# Optional: gene ID→symbol map (two columns: gene_id  gene_name).
# If you don't have one, provide a GTF and we'll build it:
ID_SYMBOL_MAP="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_1Mbp_window/00_prepare_genome/03_id_symbol_map/id_symbol_map.txt"     # or leave empty
GTF="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf"     # used only if ID_SYMBOL_MAP=""

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
### ========================

mkdir -p ${WORK_DIR}

# 12)
python3 ${PROJ_DIR}/12_plot_global_tf_ranking.py \
  --tf_rank_annot "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.tf_ranking.annot.tsv" \
  --coefs         "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.filtered.tsv" \
  --topN 20 \
  --flip-coef \
  --out "${WORK_DIR}/fig_global_tf_ranking_top20.pdf"

# 13) 
python3 ${PROJ_DIR}/13_plot_tf_gene_heatmap.py \
  --summary        "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.topTF_topGenes.summary.with_gabc.tsv" \
  --tf_rank_annot  "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.tf_ranking.annot.tsv" \
  --topT 3 --topG 25 \
  --out "${WORK_DIR}/fig_tf_gene_heatmap_score_T10xG20.pdf"

# 14)
mkdir -p "${WORK_DIR}/tf_panels"
python3 ${PROJ_DIR}/14_plot_per_tf_panels.py \
  --summary        "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.topTF_topGenes.summary.with_gabc.tsv" \
  --tf_rank_annot  "${WORK_DIR}/Old_LTHSC:Young_LTHSC_ATAC-seq.tf_ranking.annot.tsv" \
  --topTF 3 --topG 25 \
  --flip-lfc --flip-coef \
  --label-top 10 --label-metric score \
  --outdir "${WORK_DIR}/tf_panels"

