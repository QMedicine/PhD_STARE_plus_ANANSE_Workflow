# Make a proper BED6 (name, score, strand filled)
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR,0,"."}' \
  /vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed/union_peaks.bed \
  > /vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed/union_peaks_6_cols.bed
