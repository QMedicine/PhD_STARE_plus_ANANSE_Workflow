#!/bin/bash
#SBATCH --job-name=Add_activity_columns
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/03_Add_activity_columns_to_union_peaks_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/03_Add_activity_columns_to_union_peaks_error=%x_%j.err

export LC_ALL=C

# ---------- Required tools ----------
module load bedtools/2.31.1

# Edit these to your BAM paths
Y_BAM1=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/SRR13186609_Young_HSC_B1_Rep1_filtered_no_blacklist.bam
Y_BAM2=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/SRR13186610_Young_HSC_B1_Rep2_filtered_no_blacklist.bam
O_BAM1=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/SRR13186627_Old_HSC_B1_Rep1_filtered_no_blacklist.bam
O_BAM2=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/SRR13186628_Old_HSC_B1_Rep2_filtered_no_blacklist.bam

# Paths
UNION="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed/union_peaks.named.bed"
OUT_TSV="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed/union_peaks_with_activity.tsv"

# Helper: counts per peak (peak_id, count), sorted by peak_id using version sort
counts() {
  local bam="$1" out="$2"
  bedtools coverage -a "$UNION" -b "$bam" -counts \
    | awk 'BEGIN{OFS="\t"}{print $4,$NF}' \
    | sort -k1,1V > "$out"
}

# Make per-replicate count tables
counts "$Y_BAM1" y1.c
counts "$Y_BAM2" y2.c
counts "$O_BAM1" o1.c
counts "$O_BAM2" o2.c

# Normalize to CPM (tab-separated), keep sorting
norm_cpm() {
  local in="$1" out="$2"
  local S
  S=$(awk '{s+=$2}END{print s}' "$in")
  awk -v S="$S" 'BEGIN{OFS="\t"}{printf "%s\t%.6f\n",$1, ($2*1e6)/S}' "$in" \
    | sort -k1,1V > "$out"
}

norm_cpm y1.c y1.norm
norm_cpm y2.c y2.norm
norm_cpm o1.c o1.norm
norm_cpm o2.c o2.norm

# Average replicates per condition
join -t $'\t' -1 1 -2 1 y1.norm y2.norm \
  | awk -F'\t' 'BEGIN{OFS="\t"}{printf "%s\t%.6f\n",$1, ($2+$3)/2}' \
  | sort -k1,1V > young.avg

join -t $'\t' -1 1 -2 1 o1.norm o2.norm \
  | awk -F'\t' 'BEGIN{OFS="\t"}{printf "%s\t%.6f\n",$1, ($2+$3)/2}' \
  | sort -k1,1V > old.avg

# Join activities back onto the BED (sort BED by peak_id with version sort)
join -t $'\t' -1 4 -2 1 <(sort -t $'\t' -k4,4V "$UNION") young.avg \
 | join -t $'\t' -1 1 -2 1 - old.avg \
 | awk -F'\t' 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' \
 > "$OUT_TSV"

# ---- Add header (so STARE_ABCpp uses column names instead of c5/c6) ----
# If the file already has the exact header, do nothing (idempotent).
HDR=$'#chrom\tstart\tend\tpeak_id\tYoung_HSC\tOld_HSC'
if ! head -n1 "$OUT_TSV" | grep -qxF "$HDR"; then
  tmp="${OUT_TSV}.with_header.tmp"
  { printf "%s\n" "$HDR"; cat "$OUT_TSV"; } > "$tmp"
  mv "$tmp" "$OUT_TSV"
  echo "[OK] Added header to $OUT_TSV"
else
  echo "[OK] Header already present in $OUT_TSV"
fi

# QC: check column count = 6 everywhere
bad=$(awk -F'\t' 'NF!=6{c++}END{print c+0}' "$OUT_TSV")
if [[ "$bad" -ne 0 ]]; then
  echo "[WARN] $bad lines did not have 6 columns in $OUT_TSV" >&2
else
  echo "[OK] Built $OUT_TSV with 6 columns (Young_HSC, Old_HSC activities present)."
fi

# Clean
rm -f y*.c o*.c *.norm young.avg old.avg
