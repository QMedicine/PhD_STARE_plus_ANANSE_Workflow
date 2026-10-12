#!/bin/bash
#SBATCH --job-name=Add_peak_width_offset_CPM
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/03_Add_peak_width_offset_CPM_to_union_summit_peaks_build_bed_%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/03_Add_peak_width_offset_CPM_to_union_summit_peaks_build_bed_%x_%j.err

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
# Loess offset normalization & CPM (NEW FILE; do not overwrite TMM output)
# =========================
module load R/4.5.1
export R_LIBS_USER="/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1"

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(csaw)
})

outdir   <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed"
count_fp <- file.path(outdir, "union_counts.tsv")
base_fp  <- file.path(outdir, "union_base.tsv")

# NEW outputs (distinct names)
out_tsv  <- file.path(outdir, "union_summit500.with_activity.loessOffsets.tsv")
qcfp     <- file.path(outdir, "loessOffsets_qc.tsv")

# Read counts
C <- fread(count_fp)
stopifnot(colnames(C)[1] == "peak_id")

mat <- as.matrix(C[, -1, with=FALSE])
storage.mode(mat) <- "integer"
colnames(mat) <- colnames(C)[-1]

# Build SummarizedExperiment with totals
se <- SummarizedExperiment(assays = list(counts = mat))
colData(se)$totals <- colSums(mat)

# (Optional; not required for normalization)
colData(se)$group <- factor(c("Young","Young","Old","Old"))

# Compute loess offsets (feature-by-sample)
se <- normOffsets(se)   # loess-based by default

# Convert to CPM-like values using offsets (linear CPM, not logCPM)
cpm_mat <- calculateCPM(se, use.offsets=TRUE, use.norm.factors=FALSE, log=FALSE, assay.id="counts")
cpm_dt  <- as.data.table(cpm_mat)
setnames(cpm_dt, colnames(cpm_mat), colnames(C)[-1])

# Combine with base BED
base <- fread(base_fp, col.names = c("#chrom","start","end","peak_id"))
res  <- cbind(base, cpm_dt)

# Condition means (assumes these exact sample colnames exist)
res[, Young_mean := (Young_rep1 + Young_rep2)/2 ]
res[, Old_mean   := (Old_rep1   + Old_rep2)/2   ]

# Write out
fwrite(res, out_tsv, sep="\t", quote=FALSE)

# QC summary for offsets (not the full matrix)
off <- assay(se, "offset")
qc <- data.table(
  sample = colnames(mat),
  lib.size = colData(se)$totals,
  offset_median = apply(off, 2, median),
  offset_IQR = apply(off, 2, function(x) diff(stats::quantile(x, c(0.25, 0.75))))
)
fwrite(qc, qcfp, sep="\t", quote=FALSE)

cat("[OK] Wrote:", out_tsv, "and", qcfp, "\n")
RSCRIPT

# =========================
# Compare normalization outputs (plots)
# =========================
module load R/4.5.1
export R_LIBS_USER="/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1"

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

outdir <- "/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed"

tmm_fp   <- file.path(outdir, "union_summit500.with_activity.tsv")  # existing (DO NOT REWRITE)
loess_fp <- file.path(outdir, "union_summit500.with_activity.loessOffsets.tsv")

plotdir <- file.path(outdir, "normalization_comparison_TMM_vs_loessOffsets")
dir.create(plotdir, showWarnings=FALSE, recursive=TRUE)

tmm   <- fread(tmm_fp)
loess <- fread(loess_fp)

# Basic checks
stopifnot(nrow(tmm) == nrow(loess))
stopifnot(all(tmm$peak_id == loess$peak_id))

# Columns (assumes these exist)
sample_cols <- c("Young_rep1","Young_rep2","Old_rep1","Old_rep2")
stopifnot(all(sample_cols %in% names(tmm)))
stopifnot(all(sample_cols %in% names(loess)))
stopifnot(all(c("Young_mean","Old_mean") %in% names(tmm)))
stopifnot(all(c("Young_mean","Old_mean") %in% names(loess)))

# Helper: compute mean-mean MA and logFC with a small prior to avoid infinities
prior <- 0.5

make_ma <- function(dt, label){
  data.table(
    method = label,
    peak_id = dt$peak_id,
    A = 0.5 * (log2(dt$Young_mean + prior) + log2(dt$Old_mean + prior)),
    M = log2(dt$Old_mean + prior) - log2(dt$Young_mean + prior)
  )
}

ma_tmm   <- make_ma(tmm,   "TMM_CPM")
ma_loess <- make_ma(loess, "LoessOffset_CPM")
ma_all   <- rbind(ma_tmm, ma_loess)

# 1) MA plot (group means)
p_ma <- ggplot(ma_all, aes(x=A, y=M)) +
  geom_point(alpha=0.15, size=0.4) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~method, ncol=2) +
  labs(title="MA plot (Young_mean vs Old_mean) by normalization",
       x="A = mean log2(CPM+prior)",
       y="M = log2FC (Old/Young)")

ggsave(file.path(plotdir, "MAplot_groupMeans_TMM_vs_LoessOffsets.pdf"),
       p_ma, width=11, height=5)

# 2) Density of log2FC (group means)
p_den <- ggplot(ma_all, aes(x=M)) +
  geom_density() +
  facet_wrap(~method, ncol=1) +
  labs(title="Distribution of log2FC (Old/Young) by normalization",
       x="log2FC (Old/Young) from mean CPM", y="Density")

ggsave(file.path(plotdir, "log2FC_density_TMM_vs_LoessOffsets.pdf"),
       p_den, width=8, height=6)

# 3) Difference-in-logFC (Loess - TMM) vs abundance (to show systematic shifts)
delta <- data.table(
  peak_id = tmm$peak_id,
  A = ma_tmm$A,  # abundance axis (either is fine; they’re very similar)
  dM = ma_loess$M - ma_tmm$M
)

p_dm <- ggplot(delta, aes(x=A, y=dM)) +
  geom_point(alpha=0.15, size=0.4) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(title="Δlog2FC = (LoessOffset - TMM) vs abundance",
       x="A = mean log2(CPM+prior)",
       y="Δlog2FC")

ggsave(file.path(plotdir, "DeltaLog2FC_LoessMinusTMM_vs_Abundance.pdf"),
       p_dm, width=8, height=5)

# 4) Per-sample concordance scatter: log10(CPM+1) TMM vs Loess
mk_scatter <- function(samp){
  dt <- data.table(
    sample = samp,
    x = log10(tmm[[samp]] + 1),
    y = log10(loess[[samp]] + 1)
  )
  r <- suppressWarnings(cor(dt$x, dt$y, method="spearman"))
  ggplot(dt, aes(x=x, y=y)) +
    geom_point(alpha=0.10, size=0.35) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    labs(title=paste0("TMM vs LoessOffset: ", samp, " (Spearman ρ=", signif(r, 3), ")"),
         x="log10(TMM_CPM + 1)",
         y="log10(LoessOffset_CPM + 1)")
}

for (s in sample_cols){
  ggsave(file.path(plotdir, paste0("scatter_log10CPM_TMM_vs_LoessOffsets_", s, ".pdf")),
         mk_scatter(s), width=6, height=5)
}

# Write a small summary table too
summ_fp <- file.path(plotdir, "summary_metrics.tsv")
summary_dt <- data.table(
  metric = c("spearman_log2FC_groupMeans", "median_abs_delta_log2FC"),
  value  = c(
    suppressWarnings(cor(ma_tmm$M, ma_loess$M, method="spearman")),
    median(abs(delta$dM), na.rm=TRUE)
  )
)
fwrite(summary_dt, summ_fp, sep="\t")

cat("[OK] Plots written to:", plotdir, "\n")
cat("[OK] Summary metrics:", summ_fp, "\n")
RSCRIPT
