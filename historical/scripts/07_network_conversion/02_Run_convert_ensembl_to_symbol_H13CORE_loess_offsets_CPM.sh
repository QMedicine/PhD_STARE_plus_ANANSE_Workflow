#!/bin/bash
#SBATCH --job-name=02_STARE_to_Binding
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/02_Run_convert_ensembl_to_symbol_H13_loess_offsets_CPM=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/02_Run_convert_ensembl_to_symbol_H13_loess_offsets_CPM=%x_%j.err


#####################################################
####### --- USER DEFINED VARIABLES --- ##############

# -------------------- Parameters -------------------
RES=100kb                  # label only (e.g., 30kb / 50kb / 100kb)
BIN=100000                 # numeric Hi-C bin size
W_VALUE=5000000           # -w gene window (total span)
M_VALUE=5000000           # -m enhancer window

# Paths
GTF="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf"
SCRIPT="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/bin/02_convert_ensembl_to_symbol.py"

# OUT_BASE constructed to mirror your earlier convention
OUT_BASE="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM"

# Files to convert (just .tsv)
YOUNG_TSV="${OUT_BASE}/Young_HSC_STARE_to_ananse_network.tsv"
OLD_TSV="${OUT_BASE}/Old_HSC_STARE_to_ananse_network.tsv"

# Output suffix for converted files
OUT_SUFFIX="_SYMBOLS"

# Dash between TF and gene (em dash U+2014)
SEP="—"

# --- Activate conda environment ---
module load miniconda3

# Source the main conda script to initialize conda in the job's shell
#    This defines the 'conda activate' function.
#    The path is taken directly from your .bashrc file.
source /stornext/System/data/software/rhel/9/base/tools/miniconda3/latest/etc/profile.d/conda.sh

# Now, activate your specific conda environment
conda activate ananse


# ----------------------------------------------

set -euo pipefail

echo "[INFO] OUT_BASE: ${OUT_BASE}"
echo "[INFO] Using GTF: ${GTF}"
echo "[INFO] Files:"
printf '  - %s\n' "$YOUNG_TSV" "$OLD_TSV"

python --version

python "$SCRIPT" \
  --gtf "$GTF" \
  --sep "$SEP" \
  --out-suffix "$OUT_SUFFIX" \
  "$YOUNG_TSV" "$OLD_TSV"

echo "[DONE] Converted files (with ${OUT_SUFFIX} appended) are in: ${OUT_BASE}"