#!/bin/bash
#SBATCH --job-name=stare_contacts
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/01_stare_contacts_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/01_stare_contacts_error=%x_%j.err

#set -euo pipefail

# ======= USER PATHS =======
IN_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/h5_and_cooler_files_BL_removed_ICE_normalized_using_HiCExplorer"
BASE_OUT="/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_files"
CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"

echo "[*] IN_DIR   : $IN_DIR"
echo "[*] BASE_OUT : $BASE_OUT"
echo "[*] CHRS     : $CHRS"

module load miniconda3
conda activate hicexplorer

# ---- sanity checks ----
command -v hicConvertFormat >/dev/null 2>&1 || { echo "[ERROR] hicConvertFormat not found." >&2; exit 1; }
mkdir -p "$BASE_OUT"
TMP_BASE="${BASE_OUT}/tmp_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMP_BASE"
echo "[*] Using TMP_BASE: $TMP_BASE"

# ======= INPUT FILES (explicit) =======
# OLD (O_) LTHSC
OLD_30K="${IN_DIR}/O_LTHSC_30000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"
OLD_50K="${IN_DIR}/O_LTHSC_50000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"
OLD_100K="${IN_DIR}/O_LTHSC_100000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"

# YOUNG (Y_) LTHSC
YOUNG_30K="${IN_DIR}/Y_LTHSC_30000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"
YOUNG_50K="${IN_DIR}/Y_LTHSC_50000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"
YOUNG_100K="${IN_DIR}/Y_LTHSC_100000_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix.h5"

# Map to conditions and nice resolution labels
# Format: "condition label_kb path"
FILES=(
  "Old_HSC 30kb ${OLD_30K}"
  "Old_HSC 50kb ${OLD_50K}"
  "Old_HSC 100kb ${OLD_100K}"
  "Young_HSC 30kb ${YOUNG_30K}"
  "Young_HSC 50kb ${YOUNG_50K}"
  "Young_HSC 100kb ${YOUNG_100K}"
)

echo "[*] Listing input matrices:"
for entry in "${FILES[@]}"; do
  set -- $entry
  cond="$1"; res="$2"; path="$3"
  if [[ -s "$path" ]]; then
    echo "    - $cond  $res  -> $path"
  else
    echo "[WARN] Missing or empty: $path (skipping this one)" >&2
  fi
done

# ======= MAIN LOOP =======
for entry in "${FILES[@]}"; do
  set -- $entry
  COND="$1"       # Old_HSC or Young_HSC
  RESLBL="$2"     # 30kb / 50kb / 100kb (label only)
  IN_H5="$3"

  # Skip if input missing
  if [[ ! -s "$IN_H5" ]]; then
    echo "[WARN] Skipping $COND $RESLBL (file not found): $IN_H5" >&2
    continue
  fi

  OUT_DIR="${BASE_OUT}/${COND}/${RESLBL}"
  mkdir -p "$OUT_DIR"
  echo "[*] Processing: condition=$COND  res=$RESLBL"
  echo "[*] Output dir: $OUT_DIR"

  # One GI (whole genome) per file; we will split by chromosome afterwards.
  GI_PREFIX="${TMP_BASE}/${COND}_${RESLBL}_all"
  GI_TSV="${GI_PREFIX}.tsv"

  echo "  - Converting H5 -> GInteractions (whole genome)..."
  hicConvertFormat \
    --matrices "$IN_H5" \
    --inputFormat h5 \
    --outputFormat ginteractions \
    --outFileName "$GI_PREFIX"

  # The tool usually appends .tsv for ginteractions; fall back if not.
  if [[ ! -s "$GI_TSV" && -s "$GI_PREFIX" ]]; then GI_TSV="$GI_PREFIX"; fi
  if [[ ! -s "$GI_TSV" ]]; then
    echo "[ERROR] Expected GI TSV not found: ${GI_PREFIX}(.tsv)" >&2
    continue
  fi
  echo "    [ok] GI file: $GI_TSV  ($(wc -l < "$GI_TSV") lines)"

  # Split GI by chromosome and write STARE format: start1  start2  contact (tab), gzipped.
  for CHR in $CHRS; do
    # Handle potential 'chr' prefix differences in GI output (e.g., '1' vs 'chr1')
    ALTCHR="${CHR#chr}"
    STARE_OUT="${OUT_DIR}/${CHR}.txt.gz"

    echo "    - ${COND} ${RESLBL} ${CHR}: filtering GI -> STARE format..."
    awk -v chr="$CHR" -v alt="$ALTCHR" 'BEGIN{OFS="\t"}
         ($1==chr || $1==alt) && ($4==chr || $4==alt) {print $2,$5,$7}' "$GI_TSV" \
      | gzip -c > "$STARE_OUT"

    if [[ -s "$STARE_OUT" ]]; then
      LINES=$(zcat "$STARE_OUT" | wc -l || echo 0)
      echo "      [ok] Wrote ${STARE_OUT}  (lines: $LINES)"
    else
      echo "      [WARN] No data for ${CHR} -> ${STARE_OUT}" >&2
    fi
  done

  # Optional: remove intermediate GI to save space
  rm -f "$GI_TSV" 2>/dev/null || true
done

echo "[*] Done. Per-chromosome STARE files are under: $BASE_OUT/<Old_HSC|Young_HSC>/<30kb|50kb|100kb>/chrN.txt.gz"
echo "[*] Remember STARE -k must match the chosen resolution (30000, 50000, or 100000)."
