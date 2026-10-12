#!/bin/bash
#SBATCH --job-name=05_run_STARE_50_50k_1M
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/05_run_STARE_50kbp_HiC_50kbp_w_1Mbp_m_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/logs/05_run_STARE_50kbp_HiC_50kbp_w_1Mbp_m_error=%x_%j.err

set -euo pipefail

module load apptainer/1.3.5

# -------------------- Container --------------------
SIF=/vast/projects/Immune_Aging/apptainer_singularity_cache/stare-abc_1.0.5.sif

# -------------------- Parameters -------------------
RES=50kb                  # label only (e.g., 30kb / 50kb / 100kb)
BIN=50000                 # numeric Hi-C bin size
W_VALUE=50000           # -w gene window (total span)
M_VALUE=1000000           # -m enhancer window

# -------------------- Inputs -----------------------
BED=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed/union_summit500.with_activity.tsv
# Columns in BED (as you built them):
# 1:chr 2:start 3:end 4:peak_id 5:Young_rep1 6:Young_rep2 7:Old_rep1 8:Old_rep2 9:Young_mean 10:Old_mean

GTF=/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf
FASTA=/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.dna.primary_assembly.fa

MOTIFS_TF=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_expanded_motif_dimers/02_motifs/06_psem/motifs.psem  # use with -p

CONTACT_Y=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_files/Young_HSC/${RES}
CONTACT_O=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_files/Old_HSC/${RES}

# -------------------- Outputs ----------------------
OUT_BASE=/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_${RES}_w${W_VALUE}_m${M_VALUE}
mkdir -p "$OUT_BASE"   # parent only; STARE will create the final subdirs

OUT_Y=${OUT_BASE}/Young_LTHSC
OUT_O=${OUT_BASE}/Old_LTHSC

# -------------------- Preflight --------------------
fail=0
[[ -r "$BED"   ]] || { echo "[ERR] BED not readable:  $BED"; fail=1; }
[[ -r "$GTF"   ]] || { echo "[ERR] GTF not readable:  $GTF"; fail=1; }
[[ -r "$FASTA" ]] || { echo "[ERR] FASTA not readable:$FASTA"; fail=1; }
[[ -r "$MOTIFS_TF" ]] || { echo "[ERR] Motifs file not readable: $MOTIFS_TF"; fail=1; }
[[ -d "$CONTACT_Y" ]] || { echo "[ERR] Contact folder missing: $CONTACT_Y"; fail=1; }
[[ -d "$CONTACT_O" ]] || { echo "[ERR] Contact folder missing: $CONTACT_O"; fail=1; }

# Contact file sanity: require gz and three columns
for testdir in "$CONTACT_Y" "$CONTACT_O"; do
  tst=$(ls -1 "$testdir"/*chr*.gz 2>/dev/null | head -n1 || true)
  if [[ -z "$tst" ]]; then
    echo "[ERR] No gzipped contact files with 'chr' in name in $testdir"
    fail=1
  else
    if ! zcat "$tst" | awk -F'\t' 'NR==1{ok=(NF==3 && $1 ~ /^[0-9XYM]+$/ && $2 ~ /^[0-9XYM]+$/ && $3 ~ /^[-+0-9.eE]+$/)} END{exit(ok?0:1)}'
    then
      echo "[ERR] Contact format check failed for: $tst  (expect: bin1<TAB>bin2<TAB>contact)"
      fail=1
    fi
  fi
done
[[ $fail -eq 0 ]] || exit 1

# Preflight (host side)
[[ -r "$MOTIFS_TF" ]] || { echo "[ERR] Motifs file not readable: $MOTIFS_TF"; fail=1; }

# -------------------- Binds ------------------------

# Binds
BED_DIR=$(dirname "$BED")
GTF_DIR=$(dirname "$GTF")
FA_DIR=$(dirname "$FASTA")
MOTIFS_DIR=$(dirname "$MOTIFS_TF")

BIND="--bind $BED_DIR,$GTF_DIR,$FA_DIR,$CONTACT_Y,$CONTACT_O,$OUT_BASE,$MOTIFS_DIR"

# Optional: one broad bind instead (if allowed)
# BIND="--bind /vast/projects/Immune_Aging"

# Test visibility inside the container
apptainer exec --cleanenv $BIND "$SIF" bash -lc "ls -l '$MOTIFS_TF' >/dev/null" \
  || { echo '[ERR] Motif file not visible inside container'; exit 1; }

# Test write
apptainer exec --cleanenv $BIND "$SIF" bash -lc "echo ok > '$OUT_BASE/_write_test'"; rm -f "$OUT_BASE/_write_test"

# -------------------- Helper: one STARE run --------------------
run_stare () {
  local CONTACT_DIR="$1"
  local N_SPEC="$2"         # e.g., "5-6" or "7-8" or "9" etc.
  local OUT_PREFIX="$3"     # e.g., $OUT_Y or $OUT_O

  echo "[*] Running STARE: contacts=$CONTACT_DIR  -n $N_SPEC  -o $OUT_PREFIX"

  # IMPORTANT: Do NOT pre-create $OUT_PREFIX; STARE will create it and abort if it exists.
  apptainer exec --cleanenv $BIND "$SIF" bash -lc "
    set -euo pipefail
    OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
    STARE.sh \
      -b '$BED' \
      -n '$N_SPEC' \
      -a '$GTF' \
      -g '$FASTA' \
      -f '$CONTACT_DIR' \
      -k $BIN \
      -w $W_VALUE \
      -m $M_VALUE \
      -t 0.02 \
      -p '$MOTIFS_TF' \
      -c ${SLURM_CPUS_PER_TASK:-32} \
      -i all_tss \
      -q True \
      -d True \
      -o '$OUT_PREFIX'
  "
}

# -------------------- Launch both conditions -------------------
# Young = columns 5-6; Old = columns 7-8. (You can also run condition means: Young=9, Old=10)
run_stare "$CONTACT_Y" "5,6" "$OUT_Y"
run_stare "$CONTACT_O" "7,8" "$OUT_O"

echo "[✓] Done: $OUT_BASE"
