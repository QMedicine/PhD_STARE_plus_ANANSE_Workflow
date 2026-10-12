#!/usr/bin/env bash
#SBATCH --job-name=phd_history_audit_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5 samtools bedtools/2.31.1
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/audit_v01
cd "$TASK_ROOT"
printf 'JOB=%s HOST=%s ARCH=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$(uname -m)"
apptainer --version; samtools --version; bedtools --version
STARE_SIF=/vast/projects/Immune_Aging/apptainer_singularity_cache/stare-abc_1.0.5.sif
sha256sum "$STARE_SIF" > stare_sif.sha256
apptainer inspect --json "$STARE_SIF" > stare_sif.inspect.json
apptainer exec --cleanenv "$STARE_SIF" bash -c 'command -v STARE.sh; command -v STARE; ls /usr/local/bin' > stare_container_inventory.txt
WINDOWS=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed/union_summit500.named.bed
head -n 3 "$WINDOWS" > first_windows.bed
for sample in SRR13186609_Young_HSC_B1_Rep1 SRR13186610_Young_HSC_B1_Rep2 SRR13186627_Old_HSC_B1_Rep1 SRR13186628_Old_HSC_B1_Rep2; do
 bam=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/${sample}_filtered_no_blacklist.bam
 samtools view -H "$bam" > "${sample}.header.sam"
 samtools view -b "$bam" '1:3051000-3053000' > "${sample}.region.bam"
 bedtools coverage -a first_windows.bed -b "${sample}.region.bam" -counts > "${sample}.counts.tsv"
done
for bam in /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/STAR_aligned_bams/*HSC*Aligned.sortedByCoord.out.bam; do
 samtools view -H "$bam" > "$(basename "$bam").header.sam"
done
printf 'AUDIT_COMPLETE\n'
