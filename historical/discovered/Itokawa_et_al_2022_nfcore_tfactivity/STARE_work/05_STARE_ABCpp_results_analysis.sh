#!/bin/bash
#SBATCH --job-name=05_STARE_ABCpp_results_analysis
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/05_STARE_ABCpp_results_analysis_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/05_STARE_ABCpp_results_analysis_error=%x_%j.err

export LC_ALL=C

###################################################################
# -- Part A) Pair-level ΔABC with drivers (Δactivity, Δcontact) --#
###################################################################

ABCpp_dir=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/ABC_outputs_50kb

In_Out_Dir=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/ABC_outputs_50kb/Analysis_of_STARE_ABC_results
mkdir -p "$In_Out_Dir"

Y="$ABCpp_dir"/Young_HSC_ABCpp_scoredInteractions_Young_HSC.txt.gz
O="$ABCpp_dir"/Old_HSC_ABCpp_scoredInteractions_Old_HSC.txt.gz

# Build comparable tables keyed by EnsemblID|PeakID
zcat "$Y" | awk -F'\t' 'NR>1{printf "%s|%s\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%.6f\n",$4,$6,$7,$8,$9,$10,$12,$13}' | sort -k1,1 > "$In_Out_Dir"/Y.pairs
# cols: key  y_signalValue y_Contact y_adaptedActivity y_scaledContact y_TSSdist y_ABC

zcat "$O" | awk -F'\t' 'NR>1{printf "%s|%s\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%.6f\n",$4,$6,$7,$8,$9,$10,$12,$13}' | sort -k1,1 > "$In_Out_Dir"/O.pairs

# Join and compute deltas
join -t $'\t' -1 1 -2 1 "$In_Out_Dir"/Y.pairs "$In_Out_Dir"/O.pairs > "$In_Out_Dir"/pairs.join
# join columns: key, y_sig, y_con, y_adapt, y_scon, y_dist, y_abc, o_sig, o_con, o_adapt, o_scon, o_dist, o_abc

awk -F'\t' 'BEGIN{OFS="\t"}{
  split($1,a,"|"); gene=a[1]; peak=a[2];
  dABC=$7-$13; dSig=$2-$8; dCon=$3-$9; dAdapt=$4-$10;
  print gene,peak,dABC,dSig,dCon,dAdapt,$5,$11,$6,$12
}' "$In_Out_Dir"/pairs.join > "$In_Out_Dir"/pairs.delta.tsv
# cols: EnsemblID  PeakID  ΔABC  Δsignal  Δcontact  ΔadaptedActivity  Young_scaledContact  Old_scaledContact  Young_TSSdist  Old_TSSdist

# Top 25 gains/losses
sort -k3,3rn "$In_Out_Dir"/pairs.delta.tsv | head -25 > "$In_Out_Dir"/top_deltaABC_gains.tsv
sort -k3,3n "$In_Out_Dir"/pairs.delta.tsv | head -25 > "$In_Out_Dir"/top_deltaABC_losses.tsv

zcat "$Y" | awk -F'\t' 'NR>1 && !seen[$4]++{print $4"\t"$5}' | sort -k1,1 > "$In_Out_Dir"/gene.map
join -t $'\t' -1 1 -2 1 <(sort -k1,1 "$In_Out_Dir"/pairs.delta.tsv) "$In_Out_Dir"/gene.map | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$11,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > "$In_Out_Dir"/pairs.delta.withGene.tsv
# EnsemblID  GeneName  PeakID  ΔABC  Δsignal  Δcontact  ΔadaptedActivity  Young_scaledC  Old_scaledC  Young_dist  Old_dist

###################################################################
#-- Part B) Gene-level changes (from GeneInfo) --#
###################################################################

GY="$ABCpp_dir"/Young_HSC_GeneInfo_Young_HSC.txt.gz
GO="$ABCpp_dir"/Old_HSC_GeneInfo_Old_HSC.txt.gz

zcat "$GY" | awk -F'\t' 'NR>1{printf "%s\t%s\t%s\t%.6f\t%.6f\t%.1f\t%d\n",$1,$2,$3,$6,$7,$8,$5}' | sort -k1,1 > "$In_Out_Dir"/Y.gene
# EnsemblID  GeneName  chr  AvgAct  AvgContact  AvgDistance  Nenh

zcat "$GO" | awk -F'\t' 'NR>1{printf "%s\t%s\t%s\t%.6f\t%.6f\t%.1f\t%d\n",$1,$2,$3,$6,$7,$8,$5}' | sort -k1,1 > "$In_Out_Dir"/O.gene

join -t $'\t' -1 1 -2 1 "$In_Out_Dir"/Y.gene "$In_Out_Dir"/O.gene > "$In_Out_Dir"/gene.join

awk -F'\t' 'BEGIN{OFS="\t"}{
  id=$1; nameY=$2; chrY=$3;
  dN=$7-$14; dAct=$4-$11; dCon=$5-$12; dDist=$6-$13;
  print id,nameY,chrY,dN,dAct,dCon,dDist
}' "$In_Out_Dir"/gene.join | sort -k5,5rn > "$In_Out_Dir"/gene.delta.tsv
# cols: EnsemblID  GeneName  chr  Δ#Enhancer  ΔAvgActivity  ΔAvgContact  ΔAvgDistance

###################################################################
#-- Part C) R analysis --#
###################################################################

module load R/4.5.1
Rscript -e 'd<-read.table("'"$In_Out_Dir"'/pairs.delta.tsv",sep="\t"); colnames(d)<-c("gene","peak","dABC","dSig","dCon","dAdapt","ySC","oSC","yDist","oDist"); cat("cor(dABC,dSig) =",cor(d$dABC,d$dSig,use="complete.obs"),"\n"); cat("cor(dABC,dCon) =",cor(d$dABC,d$dCon,use="complete.obs"),"\n");' > "$In_Out_Dir"/correlation_results.txt