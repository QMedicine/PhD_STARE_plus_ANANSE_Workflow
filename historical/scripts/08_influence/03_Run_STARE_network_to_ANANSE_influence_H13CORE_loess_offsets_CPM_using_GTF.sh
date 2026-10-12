#!/bin/bash
#SBATCH --job-name=03_STARE_network_to_influence
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/03_Run_STARE_network_to_ANANSE_influence_H13CORE_loess_offsets_CPM_using_GTF_=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/logs/03_Run_STARE_network_to_ANANSE_influence_H13CORE_loess_offsets_CPM_using_GTF_=%x_%j.err

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
OUT_BASE="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_${RES}_w${W_VALUE}_m${M_VALUE}_H13_motifsV2_loess_offsets_CPM"
mkdir -p "${OUT_BASE}"

OUTPUT_DIR_O2Y="${OUT_BASE}/TF_Influence_Old2Young"
mkdir -p "${OUTPUT_DIR_O2Y}"
OUTPUT_DIR_Y2O="${OUT_BASE}/TF_Influence_Young2Old"
mkdir -p "${OUTPUT_DIR_Y2O}"

# Remember tfactivity comparison is always B/A
DE_GENES_Y_over_O="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv"

DE_GENES_O_over_Y="/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv"

REF_GENOME="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa"

GTF_FILE="/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf"

#####################################
# --- Activate conda environment ---#

module load miniconda3

# Source the main conda script to initialize conda in the job's shell
#    This defines the 'conda activate' function.
#    The path is taken directly from your .bashrc file.
source /stornext/System/data/software/rhel/9/base/tools/miniconda3/latest/etc/profile.d/conda.sh

# Now, activate your specific conda environment
conda activate ananse

#############################################
# --- GimmeMotifs Cache Workaround START ---#
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


#####################################################
########## --- ANANSE SCRIPT --- ####################

# remove the -a ${GTF_FILE} \ argument as the gene symbols should match

# Target -> Young
ananse influence  -s ${OUT_BASE}/Old_HSC_STARE_to_ananse_network_SYMBOLS.tsv \
                    -t ${OUT_BASE}/Young_HSC_STARE_to_ananse_network_SYMBOLS.tsv \
                    -d ${DE_GENES_Y_over_O} \
                    -a ${GTF_FILE} \
                    -i 100000 \
                    -f \
                    -o ${OUTPUT_DIR_O2Y}/STARE_OldHSC_to_YoungHSC_TF_influence.txt \
                    -n 16

echo "finished scoring Target -> Young ananse TF influence"

# Target -> Old
ananse influence  -s ${OUT_BASE}/Young_HSC_STARE_to_ananse_network_SYMBOLS.tsv \
                    -t ${OUT_BASE}/Old_HSC_STARE_to_ananse_network_SYMBOLS.tsv \
                    -d ${DE_GENES_O_over_Y} \
                    -a ${GTF_FILE} \
                    -i 100000 \
                    -f \
                    -o "${OUTPUT_DIR_Y2O}/STARE_YoungHSC_to_OldHSC_TF_influence.txt" \
                    -n 16

echo "finished scoring Target -> Old ananse TF influence"

# --- Cleanup ---
# Good practice to remove the temporary cache directory and its contents.
echo "Cleaning up temporary cache: ${NEW_CACHE}"
rm -rf ${NEW_CACHE}

conda deactivate