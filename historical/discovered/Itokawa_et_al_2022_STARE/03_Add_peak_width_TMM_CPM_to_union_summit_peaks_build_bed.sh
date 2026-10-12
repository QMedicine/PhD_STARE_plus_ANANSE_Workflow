#!/bin/bash
#SBATCH --job-name=Add_peak_width_TMM_CPM
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/03_Add_peak_width_TMM_CPM_to_union_summit_peaks_build_bed_%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/03_Add_peak_width_TMM_CPM_to_union_summit_peaks_build_bed_%x_%j.err

set -euo pipefail

# =========================
# Paths & inputs
# =========================
BAMDIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered
OUTDIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed
mkdir -p "$OUTDIR"

# Your fixed-width (500 bp), non-merged windows (Option A)
WINDOWS_500="$OUTDIR/all_summits.500bp.clean.bed"   # 3-column BED (no IDs yet)

# Named windows: chr start end peak_id (keeps input order)
BEDWIN="$OUTDIR/union_summit500.named.bed"

# Replicate BAMs (order => columns)
declare -a BAMS=(
  SRR13186609_Young_HSC_B1_Rep1_filtered_no_blacklist.bam
  SRR13186610_Young_HSC_B1_Rep2_filtered_no_blacklist.bam
  SRR13186627_Old_HSC_B1_Rep1_filtered_no_blacklist.bam
  SRR13186628_Old_HSC_B1_Rep2_filtered_no_blacklist.bam
)

# Column labels (same order as BAMS)
declare -a LABELS=("Young_rep1" "Young_rep2" "Old_rep1" "Old_rep2")

# =========================
# Tools
# =========================
module load bedtools/2.31.1

# =========================
# Prepare named windows (no merge; exact 500 bp)
# =========================
if [[ ! -s "$BEDWIN" ]]; then
  if [[ ! -s "$WINDOWS_500" ]]; then
    echo "[ERROR] Missing 500bp windows: $WINDOWS_500" >&2
    exit 1
  fi
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "$WINDOWS_500" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR}' \
    > "$BEDWIN"
  echo "[OK] Built named windows: $BEDWIN  (n=$(wc -l < "$BEDWIN"))"
else
  echo "[OK] Using existing named windows: $BEDWIN  (n=$(wc -l < "$BEDWIN"))"
fi

# Base columns for the final table
BASE="$OUTDIR/union_base.tsv"
cut -f1-4 "$BEDWIN" > "$BASE"

# =========================
# Raw counts per window per BAM
# =========================
# bedtools coverage -counts reports, for each -a interval, the number of alignments from -b overlapping it.
# We rely on the order of -a (BEDWIN) being preserved -> columns remain aligned across replicates.

COUNTMAT="$OUTDIR/union_counts.tsv"
TMP_COUNTS_DIR="$OUTDIR/tmp_counts"
mkdir -p "$TMP_COUNTS_DIR"

# Collect all count files
for i in "${!BAMS[@]}"; do
  bam="${BAMS[$i]}"
  lab="${LABELS[$i]}"
  bam_path="$BAMDIR/$bam"

  [[ -s "$bam_path" ]] || { echo "[ERROR] Missing BAM: $bam_path" >&2; exit 1; }

  echo "[*] Counting overlaps for $lab ..."
  # Output: same BED as -a plus extra last column = count
  bedtools coverage -a "$BEDWIN" -b "$bam_path" -counts > "$TMP_COUNTS_DIR/${lab}.counts.bed"

  # Extract counts column only (last field) and save to a temporary file
  awk '{print $NF}' "$TMP_COUNTS_DIR/${lab}.counts.bed" > "$TMP_COUNTS_DIR/${lab}.counts"
done

# Assemble count matrix with header
{
  printf "peak_id"
  for lab in "${LABELS[@]}"; do printf "\t%s" "$lab"; done
  printf "\n"
  
  # Use paste on a list of files to assemble the data columns correctly
  paste <(cut -f4 "$BEDWIN") $(find "$TMP_COUNTS_DIR" -name "*.counts" | sort)
} > "$COUNTMAT"

# Cleanup
rm -rf "$TMP_COUNTS_DIR"

echo "[OK] Count matrix: $COUNTMAT"

# =========================
# TMM normalisation & CPM in R
# =========================
module load R/4.5.1
export R_LIBS_USER="/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1"

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

outdir <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed"
count_fp <- file.path(outdir, "union_counts.tsv")
base_fp  <- file.path(outdir, "union_base.tsv")
out_tsv  <- file.path(outdir, "union_summit500.with_activity.tsv")

# Read counts
C <- fread(count_fp)
stopifnot(colnames(C)[1] == "peak_id")
pid <- C[[1]]
mat <- as.matrix(C[ , -1, with=FALSE])
storage.mode(mat) <- "integer"

# edgeR DGEList
grp <- factor(c("Young","Young","Old","Old"))
dge <- DGEList(counts=mat, group=grp)

# Filter extremely low windows? (Optional; we keep all to match STARE union)
# dge <- dge[rowSums(dge$counts) > 0, ]

# TMM normalization factors
dge <- calcNormFactors(dge, method="TMM")

# CPM with normalized library sizes
cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE)  # linear scale CPM, not logCPM
cpm_dt  <- as.data.table(cpm_mat)
setnames(cpm_dt, colnames(cpm_mat), colnames(C)[-1])

# Combine with base BED
base <- fread(base_fp, col.names = c("#chrom","start","end","peak_id"))
res  <- cbind(base, cpm_dt)
# Condition means
res[, Young_mean := (Young_rep1 + Young_rep2)/2 ]
res[, Old_mean   := (Old_rep1   + Old_rep2)/2   ]

# Write out
fwrite(res, out_tsv, sep="\t", quote=FALSE)

# Also write a small QC about factors & library sizes
qcfp <- file.path(outdir, "TMM_normalization_factors.tsv")
nf <- data.table(sample=colnames(dge$counts),
                 lib.size=dge$samples$lib.size,
                 norm.factors=dge$samples$norm.factors,
                 eff.lib.size=dge$samples$lib.size * dge$samples$norm.factors)
fwrite(nf, qcfp, sep="\t")
cat("[OK] Wrote:", out_tsv, "and", qcfp, "\n")
RSCRIPT

# =========================
# QC and header check
# =========================
# Confirm expected columns
HDR=$'#chrom\tstart\tend\tpeak_id\t'"${LABELS[0]}"$'\t'"${LABELS[1]}"$'\t'"${LABELS[2]}"$'\t'"${LABELS[3]}"$'\tYoung_mean\tOld_mean'
# Ensure header exactly matches (add if needed)
if ! head -n1 "$OUTDIR/union_summit500.with_activity.tsv" | grep -qxF "$HDR"; then
  tmpH="$OUTDIR/.with_activity.header.tmp"
  { printf "%s\n" "$HDR"; tail -n +2 "$OUTDIR/union_summit500.with_activity.tsv"; } > "$tmpH"
  mv "$tmpH" "$OUTDIR/union_summit500.with_activity.tsv"
fi

# Row count consistency
n_win=$(wc -l < "$BEDWIN")
n_rows=$(($(wc -l < "$OUTDIR/union_summit500.with_activity.tsv") - 1))
if [[ "$n_win" -ne "$n_rows" ]]; then
  echo "[WARN] Row mismatch: windows=$n_win vs table=$n_rows" >&2
else
  echo "[OK] union_summit500.with_activity.tsv has $n_rows rows and expected 10 columns."
fi

echo "[✓] Done. Output: $OUTDIR/union_summit500.with_activity.tsv"
