#!/bin/bash
#SBATCH --job-name=01_STARE_to_Binding
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/01_Run_STARE_TF_Gene_Affinities_to_ananse_network_H13CORE_loess_offsets=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/01_Run_STARE_TF_Gene_Affinities_to_ananse_network_H13CORE_loess_offsets=%x_%j.err

set -euo pipefail

# Pipeline from https://github.com/vanheeringen-lab/ANANSE
# and paper:
# ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination
# Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Gert Jan C Veenstra, Huiqing Zhou, Simon J van Heeringen 
# Nucleic Acids Research, gkab598, https://doi.org/10.1093/nar/gkab598

#####################################################
####### --- USER DEFINED VARIABLES --- ##############

# -------------------- Parameters -------------------
RES=100kb                  # label only (e.g., 30kb / 50kb / 100kb)
BIN=100000                 # numeric Hi-C bin size
W_VALUE=5000000           # -w gene window (total span)
M_VALUE=5000000           # -m enhancer window

# -------------------- Outputs ----------------------
OUT_BASE=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM
mkdir -p "$OUT_BASE"

# --- STARE TF-gene affinity inputs --- #
STARE_INPUT_OLD_REP1=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM/Old_LTHSC/Gene_TF_matrices/Old_LTHSC_TF_Gene_Affinities_Old_rep1.txt.gz
STARE_INPUT_OLD_REP2=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM/Old_LTHSC/Gene_TF_matrices/Old_LTHSC_TF_Gene_Affinities_Old_rep2.txt.gz

STARE_INPUT_YOUNG_REP1=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM/Young_LTHSC/Gene_TF_matrices/Young_LTHSC_TF_Gene_Affinities_Young_rep1.txt.gz
STARE_INPUT_YOUNG_REP2=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM/Young_LTHSC/Gene_TF_matrices/Young_LTHSC_TF_Gene_Affinities_Young_rep2.txt.gz

# --- Gene expression TPM matrix --- #
EXPR_TPM="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_H13CORE_motifs/01_counts/02_tpm/counts.tpm.tsv"

IDMAP="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/00_prepare_genome/03_id_symbol_map/id_symbol_map.txt"

# --- Activate conda environment ---
module load miniconda3

# Source the main conda script to initialize conda in the job's shell
#    This defines the 'conda activate' function.
#    The path is taken directly from your .bashrc file.
source /stornext/System/data/software/rhel/9/base/tools/miniconda3/latest/etc/profile.d/conda.sh

# Now, activate your specific conda environment
conda activate ananse

# --- GimmeMotifs Cache Workaround START ---
# We need a unique cache directory for this specific job to prevent concurrent access errors.

# 1. Define a NEW_CACHE path on the fast /vast/scratch storage.
#    Using $SLURM_JOB_ID makes this directory unique for every job run.
#    (Using the path structure from your .bashrc for apptainer)
NEW_CACHE="/vast/scratch/users/$USER/gimmemotifs_cache/${SLURM_JOB_ID}"

# 2. Create this job-specific directory.
mkdir -p ${NEW_CACHE}

# 3. Export XDG_CACHE_HOME to point to this new directory.
#    GimmeMotifs will now create a new, clean cache here.
export XDG_CACHE_HOME=${NEW_CACHE}

echo "Redirected GimmeMotifs cache to job-specific directory: ${XDG_CACHE_HOME}"

# NOTE: We have REMOVED the 'cp' command. We do NOT want to copy
# the old, corrupted cache. We want GimmeMotifs to start fresh.

# --- GimmeMotifs Cache Workaround END ---

conda activate /home/users/allstaff/stacpoole.q/.conda/envs/ananse

PYTHON_SCRIPTS="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/bin"

# Young
python "${PYTHON_SCRIPTS}/01_STARE_TF_Gene_Affinities_to_ananse_network.py" \
  --stare "${STARE_INPUT_YOUNG_REP1}" "${STARE_INPUT_YOUNG_REP2}" \
  --tpm "$EXPR_TPM" \
  --outdir "${OUT_BASE}" \
  --idmap "${IDMAP}" \
  --basename Young_HSC_STARE_to_ananse_network

# Old
python "${PYTHON_SCRIPTS}/01_STARE_TF_Gene_Affinities_to_ananse_network.py" \
  --stare "${STARE_INPUT_OLD_REP1}" "${STARE_INPUT_OLD_REP2}" \
  --tpm "$EXPR_TPM" \
  --outdir "${OUT_BASE}" \
  --idmap "${IDMAP}" \
  --basename Old_HSC_STARE_to_ananse_network