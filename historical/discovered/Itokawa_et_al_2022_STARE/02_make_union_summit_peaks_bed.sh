#!/bin/bash
#SBATCH --job-name=union_peaks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/02_make_union_summit_peaks_bed_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/02_make_union_summit_peaks_bed_error=%x_%j.err

# modules you likely already have
module load bedtools/2.31.1

SUMDIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/macs3_peaks
OUTDIR=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed
mkdir -p "$OUTDIR"

# mm39 chrom sizes and blacklist/excluderanges if you have them
CHROMSIZES=/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39_ensembl.chrom.sizes
BLACKLIST=/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed  # or leave empty

# 1) Collect summits from all samples
cat \
  "$SUMDIR"/SRR13186609_Young_HSC_B1_Rep1_summits.bed \
  "$SUMDIR"/SRR13186610_Young_HSC_B1_Rep2_summits.bed \
  "$SUMDIR"/SRR13186627_Old_HSC_B1_Rep1_summits.bed \
  "$SUMDIR"/SRR13186628_Old_HSC_B1_Rep2_summits.bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1}' \
| sort -k1,1 -k2,2n \
> "$OUTDIR"/all_summits.1bp.bed

# 2) Make fixed-width windows (±250 bp around summit), clipped to genome
bedtools slop -i "$OUTDIR"/all_summits.1bp.bed -g "$CHROMSIZES" -b 250 \
| sort -k1,1 -k2,2n \
> "$OUTDIR"/all_summits.500bp.bed

# 3) Optional: remove blacklisted loci
if [[ -s "$BLACKLIST" ]]; then
  bedtools subtract -A -a "$OUTDIR"/all_summits.500bp.bed -b "$BLACKLIST" \
  > "$OUTDIR"/all_summits.500bp.clean.bed
else
  cp "$OUTDIR"/all_summits.500bp.bed "$OUTDIR"/all_summits.500bp.clean.bed
fi

