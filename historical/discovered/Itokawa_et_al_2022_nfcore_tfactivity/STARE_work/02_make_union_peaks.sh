#!/bin/bash
#SBATCH --job-name=union_peaks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/02_make_union_peaks_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/02_make_union_peaks_error=%x_%j.err

# ---------- User inputs ----------
BED_OLD_R1="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/macs3_peaks/narrowPeak_as_bed/SRR13186627_Old_HSC_B1_Rep1_peaks.bed"
BED_OLD_R2="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/macs3_peaks/narrowPeak_as_bed/SRR13186628_Old_HSC_B1_Rep2_peaks.bed"
BED_YNG_R1="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/macs3_peaks/narrowPeak_as_bed/SRR13186609_Young_HSC_B1_Rep1_peaks.bed"
BED_YNG_R2="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/macs3_peaks/narrowPeak_as_bed/SRR13186610_Young_HSC_B1_Rep2_peaks.bed"

# Optional: mm39 blacklist/excluderanges (leave empty "" to skip)
BLACKLIST="/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed"

OUTDIR="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed"
mkdir -p "$OUTDIR"/tmp

echo "[*] Output directory: $OUTDIR"

# ---------- Required tools ----------
module load bedtools/2.31.1 

# ---------- Helper: clean to BED3+sorted (and optional blacklist) ----------
clean() {
  local in="$1"; local out="$2"
  if [[ ! -s "$in" ]]; then echo "[ERROR] missing: $in"; exit 1; fi
  # Keep first 3 cols only and sort; many 'narrowPeak_as_bed' already have BED3
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$in" \
    | sort -k1,1 -k2,2n > "${out}.tmp"

  if [[ -n "${BLACKLIST}" && -s "${BLACKLIST}" ]]; then
    bedtools subtract -A -a "${out}.tmp" -b "${BLACKLIST}" > "$out"
    rm -f "${out}.tmp"
  else
    mv "${out}.tmp" "$out"
  fi
  echo "    [ok] cleaned: $in -> $out  (n=$(wc -l < "$out"))"
}

# ---------- Clean each BED ----------
C1="$OUTDIR/tmp/old_r1.clean.bed";  clean "$BED_OLD_R1" "$C1"
C2="$OUTDIR/tmp/old_r2.clean.bed";  clean "$BED_OLD_R2" "$C2"
C3="$OUTDIR/tmp/yng_r1.clean.bed";  clean "$BED_YNG_R1" "$C3"
C4="$OUTDIR/tmp/yng_r2.clean.bed";  clean "$BED_YNG_R2" "$C4"

# ---------- Union merge ----------
# bedtools merge collapses overlapping or directly adjacent ("book-ended") intervals into one.
# This is ideal for defining a non-redundant candidate set for ABC/STARE.
cat "$C1" "$C2" "$C3" "$C4" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$OUTDIR/union_peaks.bed"

echo "[*] union_peaks.bed written: $(wc -l < "$OUTDIR/union_peaks.bed") intervals"

# Optional: if you prefer to NOT merge directly adjacent peaks, use:
# ... | bedtools merge -i - -d -1 > "$OUTDIR/union_peaks.bed"

# ---------- Add peak IDs ----------
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR}' "$OUTDIR/union_peaks.bed" \
  > "$OUTDIR/union_peaks.named.bed"
echo "[*] union_peaks.named.bed (with peak IDs) written."

# ---------- Provenance: per-file support counts on the union ----------
# Appends 4 columns = counts of overlaps per input file, in the same order (C1..C4).
bedtools annotate -i "$OUTDIR/union_peaks.bed" -files "$C1" "$C2" "$C3" "$C4" -counts \
  > "$OUTDIR/union_peaks.annot.tsv"

# Add a simple header with recognizable sample names
bn1=$(basename "$BED_OLD_R1"); bn2=$(basename "$BED_OLD_R2"); bn3=$(basename "$BED_YNG_R1"); bn4=$(basename "$BED_YNG_R2")
{
  printf "#chrom\tstart\tend\tcount_%s\tcount_%s\tcount_%s\tcount_%s\n" "$bn1" "$bn2" "$bn3" "$bn4"
  cat "$OUTDIR/union_peaks.annot.tsv"
} > "$OUTDIR/union_peaks.annot.with_header.tsv"

echo "[*] union_peaks.annot.with_header.tsv written."

# ---------- Optional: per-condition support summary ----------
awk 'BEGIN{OFS="\t"} NR==1{print "#chrom","start","end","Old_support","Young_support"; next} {print $1,$2,$3,$4+$5,$6+$7}' \
  "$OUTDIR/union_peaks.annot.with_header.tsv" > "$OUTDIR/union_peaks.condition_counts.tsv"

echo "[*] union_peaks.condition_counts.tsv written."
echo "[✓] Done."
