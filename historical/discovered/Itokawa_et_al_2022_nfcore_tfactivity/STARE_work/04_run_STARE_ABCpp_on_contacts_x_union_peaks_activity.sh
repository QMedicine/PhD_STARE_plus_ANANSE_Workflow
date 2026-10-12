#!/bin/bash
#SBATCH --job-name=04_run_STARE_ABCpp
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/04_run_STARE_ABCpp_output=%x_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/logs/04_run_STARE_ABCpp_error=%x_%j.err

set -euo pipefail

# --General rationale of STARE--#
# STARE implements the Activity-By-Contact (ABC) / gABC model to link candidate
# enhancers to genes: score(e,g) ∝ Activity(e) × Contact(e,g), normalized across
# candidates around each gene’s TSS. With -q True (gABC), enhancer activity is
# also distributed across contacted genes; -i all_tss averages over all TSS.
# We run STARE once per condition to generate condition-specific ABC link tables.

module load apptainer/1.3.5

SIF=/vast/projects/Immune_Aging/apptainer_singularity_cache/stare-abc_1.0.5.sif

RES=50kb             # user-defined resolution (30kb / 50kb / 100kb)
RES_Num=50000        # numeric bin size that matches $RES

# -- Purpose of the union peaks with activity file -- #
# BED columns: 1:chrom 2:start 3:end 4:peak_id 5:Young_HSC_activity 6:Old_HSC_activity
# We select which activity to use with -n (1-based column index).
BED=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/union_bed/union_peaks_with_activity.tsv

# -- How STARE uses the GTF file -- #
# GTF provides gene models/TSS. With -i all_tss, STARE averages ABC across all
# annotated TSS per gene (recommended). Alternative: -i 5_tss for 5' TSS only.
GTF=/vast/projects/Immune_Aging/Reference_genomes/mm39/ensembl_release-113/Mus_musculus.GRCm39.113.gtf

CONTACT_Y=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/STARE_files/Young_HSC/$RES
CONTACT_O=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/STARE_files/Old_HSC/$RES

OUT_ABC=/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/ABC_outputs_$RES
mkdir -p "$OUT_ABC"

echo "[*] Using:"
echo "    BED     : $BED"
echo "    GTF     : $GTF"
echo "    CONTACTY: $CONTACT_Y"
echo "    CONTACTO: $CONTACT_O"
echo "    BIN     : $RES_Num"
echo "    OUT     : $OUT_ABC"

# ----------------- Preflight checks (added for robustness) -----------------
fail=0
[[ -r "$BED" ]] || { echo "[ERR] BED not readable: $BED"; fail=1; }
[[ -r "$GTF" ]] || { echo "[ERR] GTF not readable: $GTF"; fail=1; }
[[ -d "$CONTACT_Y" ]] || { echo "[ERR] Contact folder missing: $CONTACT_Y"; fail=1; }
[[ -d "$CONTACT_O" ]] || { echo "[ERR] Contact folder missing: $CONTACT_O"; fail=1; }
[[ -w "$OUT_ABC" ]] || { echo "[ERR] Output dir not writable: $OUT_ABC"; fail=1; }

yc=$(ls -1 "$CONTACT_Y"/*.txt* 2>/dev/null | head -n 3 | wc -l || true)
oc=$(ls -1 "$CONTACT_O"/*.txt* 2>/dev/null | head -n 3 | wc -l || true)
[[ "$yc" -ge 1 ]] || { echo "[ERR] No contact files in $CONTACT_Y"; fail=1; }
[[ "$oc" -ge 1 ]] || { echo "[ERR] No contact files in $CONTACT_O"; fail=1; }

# quick format sniff on one contact file: expect 3 cols; col1/col2 ints; col3 numeric
test_file="$(ls "$CONTACT_Y"/*.txt* | head -n1)"
if ! zcat -f "$test_file" | awk -F'\t' 'NR==1{ok=(NF==3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[-+0-9.eE]+$/)} END{exit (ok?0:1)}'; then
  echo "[ERR] Contact format check failed for $test_file (expect: bin1_start<TAB>bin2_start<TAB>contact)"; fail=1
fi

# BED should have exactly 6 columns
awk -F'\t' 'NF!=6{bad++} END{if(bad>0) exit 1}' "$BED" \
  || { echo "[ERR] BED does not have 6 columns everywhere: $BED"; fail=1; }

[[ "$fail" -eq 0 ]] || { echo "[ABORT] Fix the issues above."; exit 1; }


# --------------- Explicit container binds (HPC-safe) ---------------
BED_DIR=$(dirname "$BED")
GTF_DIR=$(dirname "$GTF")
BIND="--bind $BED_DIR,$GTF_DIR,$CONTACT_Y,$CONTACT_O,$OUT_ABC"

# write test inside container to OUT_ABC (checks bind + perms)
apptainer exec $BIND "$SIF" bash -lc "echo container_write_ok > '$OUT_ABC/_container_write_test.txt'" \
  || { echo "[ERR] Container cannot write to $OUT_ABC (bind/mount issue)."; exit 1; }
rm -f "$OUT_ABC/_container_write_test.txt"

# show what the container sees in the contact dirs (first few files)
apptainer exec $BIND "$SIF" bash -lc "echo '[inside] CONTACT_Y listing:'; ls -l '$CONTACT_Y'"
apptainer exec $BIND "$SIF" bash -lc "echo '[inside] CONTACT_O listing:'; ls -l '$CONTACT_O'"

# ---- STARE_ABCpp options (why these) ----
#  -b : candidate enhancers (with activity columns)
#  -n : which activity column to use (5=Young, 6=Old; 1-based)
#  -a : gene annotation (GTF) to locate TSS
#  -f : folder with per-chromosome normalized Hi-C files (names == chromosome)
#  -k : Hi-C bin size (must match contact files)
#  -w : ABC search window (~5 Mb total; canonical)
#  -t : write only pairs with ABC ≥ 0.02 (set 0 to output all)
#  -c : threads
#  -i : TSS mode (all_tss recommended)
#  -q : adapted ABC (gABC) – redistributes enhancer activity across genes by contact
#  -d : add pseudocount to contacts (robust to zeros/sparsity)
#  -o : output prefix
# We cd into $OUT_ABC *inside* the container and use a relative -o prefix. This avoids
# host→container path confusion and fixes “Can’t open the output path” errors.

run_stare () {
  local cond="$1" col="$2" contact="$3"
  echo "[*] Running STARE for $cond (activity column $col; contacts: $contact)"
  apptainer exec $BIND "$SIF" bash -lc "\
    set -euo pipefail; \
    cd '$OUT_ABC'; \
    OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}; \
    STARE_ABCpp \
      -b '$BED' -n $col -a '$GTF' \
      -f '$contact' -k $RES_Num -w 200000000 -t 0.02 -c ${SLURM_CPUS_PER_TASK:-32} \
      -i all_tss -q True -d True \
      -o '${cond}' \
  " || { echo '[ERR] STARE_ABCpp failed for' \"$cond\"; exit 1; }

  local outfile_can="$OUT_ABC/${cond}_ABCpp_scoredInteractions_${cond}.txt.gz"
  local outfile_col="$OUT_ABC/${cond}_ABCpp_scoredInteractions_c${col}.txt.gz"

  if [[ -s "$outfile_can" ]]; then
    echo "[OK] Wrote $outfile_can"
  elif [[ -s "$outfile_col" ]]; then
    echo "[OK] Wrote $outfile_col"
    # create a canonical alias for downstream steps
    ln -sfn "$(basename "$outfile_col")" "$outfile_can"
    echo "[OK] Symlinked to $outfile_can"
  else
    echo "[ERR] Expected output missing: $outfile_can (or fallback $outfile_col)"
    exit 1
  fi

}

# 1-based columns in BED: 5=Young activity, 6=Old activity
run_stare "Young_HSC" 5 "$CONTACT_Y"
run_stare "Old_HSC"   6 "$CONTACT_O"

# ---- Post-check: quick summaries per condition ----
summarize_abc () {
  local f="$1"
  local base=$(basename "$f")
  echo "----- Summary for $base -----"

  # Detect header (assume header if first line contains alphabetic tokens)
  local first
  first=$(zcat "$f" | head -n1 || true)
  local has_header=0
  if echo "$first" | grep -Eq '[A-Za-z]'; then has_header=1; fi

  # Total pairs (subtract header if present)
  local total
  total=$(zcat "$f" | wc -l)
  local pairs=$total
  if [[ $has_header -eq 1 ]]; then pairs=$((total-1)); fi
  echo "pairs_total: $pairs"

  # Unique enhancers (assumes first 3 columns are enhancer coords)
  if [[ $has_header -eq 1 ]]; then
    zcat "$f" | tail -n +2 | awk -F'\t' '{print $1":"$2"-"$3}' | sort -u | wc -l | awk '{print "enhancers_unique:",$1}'
  else
    zcat "$f" | awk -F'\t' '{print $1":"$2"-"$3}' | sort -u | wc -l | awk '{print "enhancers_unique:",$1}'
  fi

  # Try to detect a gene column by header name; else skip
  local gene_idx=""
  if [[ $has_header -eq 1 ]]; then
    gene_idx=$(
      echo "$first" | awk -F'\t' '{
        for (i=1;i<=NF;i++){
          l=tolower($i);
          if (l=="gene" || l=="gene_id" || l=="geneid" || l=="gene_name" || l=="symbol" || l=="targetgene" || l=="target_gene"){print i; exit}
        }
      }'
    )
  fi
  if [[ -n "$gene_idx" ]]; then
    if [[ $has_header -eq 1 ]]; then
      zcat "$f" | tail -n +2 | awk -v g="$gene_idx" -F'\t' '{print $g}' | sort -u | wc -l | awk '{print "genes_unique:",$1}'
    else
      zcat "$f" | awk -v g="$gene_idx" -F'\t' '{print $g}' | sort -u | wc -l | awk '{print "genes_unique:",$1}'
    fi
  else
    echo "genes_unique: (unknown – gene column not detected in header)"
  fi

  # ABC score stats: try to detect a score column by name; compute mean & median
  local score_idx=""
  if [[ $has_header -eq 1 ]]; then
    score_idx=$(
      echo "$first" | awk -F'\t' '{
        for (i=1;i<=NF;i++){
          l=tolower($i);
          if (l ~ /adapted.*abc/ || l=="abc" || l=="abc_score" || l=="score_abc" || l=="abcscore" || l=="score"){print i; exit}
        }
      }'
    )
  fi
  if [[ -n "$score_idx" ]]; then
    # mean
    local mean
    if [[ $has_header -eq 1 ]]; then
      mean=$(zcat "$f" | tail -n +2 | awk -v s="$score_idx" -F'\t' '{x+=$s; n++} END{if(n>0) printf "%.6f\n", x/n; else print "NA"}')
    else
      mean=$(zcat "$f" | awk -v s="$score_idx" -F'\t' '{x+=$s; n++} END{if(n>0) printf "%.6f\n", x/n; else print "NA"}')
    fi
    # median (numeric sort)
    local med
    if [[ $has_header -eq 1 ]]; then
      med=$(zcat "$f" | tail -n +2 | awk -v s="$score_idx" -F'\t' '{print $s}' | sort -g | awk '{a[NR]=$1} END{if(NR==0) print "NA"; else if (NR%2==1) printf "%.6f\n", a[(NR+1)/2]; else printf "%.6f\n", (a[NR/2]+a[NR/2+1])/2}')
    else
      med=$(zcat "$f" | awk -v s="$score_idx" -F'\t' '{print $s}' | sort -g | awk '{a[NR]=$1} END{if(NR==0) print "NA"; else if (NR%2==1) printf "%.6f\n", a[(NR+1)/2]; else printf "%.6f\n", (a[NR/2]+a[NR/2+1])/2}')
    fi
    echo "abc_mean: $mean"
    echo "abc_median: $med"
  else
    echo "abc_mean/median: (score column not detected)"
  fi

  # Peek first few data lines
  echo "[head]"
  zcat "$f" | head -n 5
}

YFILE="$OUT_ABC/Young_HSC_ABCpp_scoredInteractions_Young_HSC.txt.gz"
OFILE="$OUT_ABC/Old_HSC_ABCpp_scoredInteractions_Old_HSC.txt.gz"
summarize_abc "$YFILE"
summarize_abc "$OFILE"

# <What the output looks like and how to interpret it>
# Outputs in $OUT_ABC include, for each condition:
#   • <COND>_ABCpp_scoredInteractions_<COND>.txt.gz
#       - The main enhancer–gene table with ABC/gABC scores.
#         Typical columns include enhancer coordinates/ID, gene ID/symbol,
#         distances, activity, contact, and the final ABC score.
#       - Filtered to ABC >= 0.02 by -t (set -t 0 to get all pairs and filter later).
#   • Metadata files summarizing parameters used for the run.
# Quick checks:
#   zcat $YFILE | head
#   zcat $OFILE | head
# Expect thousands of pairs; counts scale with window size, peak density, and cutoff.
#
# Next step: can reuse these ABC links inside nf-core/tfactivity by adding a small
# Nextflow config for the STARE process, e.g.:
#   withName: /.*STARE.*/ {
#     ext.args = "-r $YFILE -w 5000000 -q True -i all_tss"
#   }
# and point BOTH conditions in your samplesheet to the SAME union BED so STARE inside
# tfactivity sees the exact same regions you used here for ABC.
echo "[✓] STARE ABC precompute + summary complete for $RES"
