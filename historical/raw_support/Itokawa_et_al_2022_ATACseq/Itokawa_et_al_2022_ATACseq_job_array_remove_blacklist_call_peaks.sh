#!/bin/bash
#SBATCH --job-name=ATAC_blacklist
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=12:00:00
#SBATCH --array=1-36
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/blacklist_rm__%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/logs/blacklist_rm_%A_%a.err

set -euo pipefail                         # fail fast and propagate pipe errors

############################
#        MODULES           #
############################
module load samtools/1.22.1
module load bedtools/2.31.1
module load MACS3/3.0.1

############################
#      USER VARIABLES      #
############################
scratch_output_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq"
project_output_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq"
blacklist_bed="/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges_ensembl_chrom_names.bed"
csv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/sample_information.csv"
peak_dir="$project_output_dir/macs3_peaks"

g_size=2654621783 #GRCm39 effective genome size

############################
#  PARSE THE SAMPLE SHEET  #
############################
header=$(head -n 1 "$csv_file")
sample_id_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Sample_ID" | cut -d: -f1)

sample_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$csv_file")
Sample_ID=$(echo "$sample_line" | awk -v c=$sample_id_col 'BEGIN{FS=","}{print $c}' | sed 's/^"\|"$//g')

############################
#   INPUT / OUTPUT FILES   #
############################
filtered_bam="$project_output_dir/filtered/${Sample_ID}_filtered.bam"
no_bl_bam="$project_output_dir/filtered/${Sample_ID}_filtered_no_blacklist.bam"
no_bl_bam_idx="${no_bl_bam}.bai"

echo "[$(date)]  Sample: $Sample_ID"
echo "Input BAM : $filtered_bam"
echo "Output BAM: $no_bl_bam"

############################
#        PIPELINE          #
############################
# 1. Remove blacklist reads and sort
bedtools intersect -v -abam "$filtered_bam" -b "$blacklist_bed" \
  | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$no_bl_bam" -

# 2. Index
samtools index -@ "${SLURM_CPUS_PER_TASK}" "$no_bl_bam" "$no_bl_bam_idx"

# 3. Determine library type
if (( $(samtools view -c -f 1 "$no_bl_bam") > 0 )); then
    libtype="paired-end"
    echo "[$(date)]  Detected paired-end BAM, will use BAMPE, etc"
    macs3 callpeak \
        -t "$no_bl_bam" \
        -f BAMPE \
        -g ${g_size} \
        -q 0.05 \
        -n "$Sample_ID" \
        --outdir "$peak_dir"
else
    libtype="single-end"
    echo "[$(date)]  Detected single-end BAM, will use nomodel, shift, extsize, etc"
    macs3 callpeak \
        -t "$no_bl_bam" \
        -g ${g_size} \
        --nomodel --shift -100 --extsize 200 \
        -q 0.05 \
        -n "$Sample_ID" \
        --outdir "$peak_dir"
fi

echo "[$(date)]  Finished $Sample_ID – $libtype processing complete"
