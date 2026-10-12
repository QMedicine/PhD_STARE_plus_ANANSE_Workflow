#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/download_sra_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/logs/download_sra_%A_%a.err
#SBATCH --array=0-27
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Load the required module
module load sra-toolkit/3.1.0

# Directory to save downloaded SRA files
download_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data/part_2_files_MPP_etc"
mkdir -p "$download_dir"

# Read the list of SRA accession IDs into an array
mapfile -t sra_ids < /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_SRR_list_p2.txt

# Use the SLURM_ARRAY_TASK_ID to select the correct accession
sra_id=${sra_ids[$SLURM_ARRAY_TASK_ID]}

echo "Downloading $sra_id..."
prefetch --max-size 2G --resume yes --output-directory "$download_dir" "$sra_id"

echo "Converting $sra_id to FASTQ..."
fasterq-dump --outdir "$download_dir" --threads 2 "$download_dir/$sra_id"

echo "$sra_id processing complete."

#Run on slurm using:
# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/download_and_convert_SRA_files_for_Itokawa_et_al_2022_RNAseq_p2.sh