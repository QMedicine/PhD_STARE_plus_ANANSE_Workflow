#!/bin/bash
#SBATCH --job-name=hic_full_valid
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/06_hic_full_valid_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/06_hic_full_valid_%A_%a.err
#SBATCH --array=0-5
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail
shopt -s nullglob

# -----------------------------------------------------------------------------
# PURPOSE
# Convert each blacklist-filtered, query-name-sorted BAM independently into a
# full-depth set of deduplicated, high-confidence Hi-C contacts.
#
# This script DOES NOT downsample or create Coolers. All six array tasks must
# finish first, because the smallest full-depth valid-pair library cannot be
# known until every sample has been processed.
# -----------------------------------------------------------------------------

BAM_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Name_sorted_merged_bam_files_BL_removed"
PROJECT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_BL_removed_coolers_pairs_files"
QC_DIR="${PROJECT_DIR}/downsampled_merged_HiC_QC"
LOG_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs"

# Full-depth valid pairs are intermediate inputs for Script 07. They are kept
# in shared scratch rather than permanent project storage.
FULL_VALID_DIR="/vast/scratch/users/stacpoole.q/HiC_full_depth_valid_pairs"
TASK_TMP="/vast/scratch/users/stacpoole.q/hic_full_valid_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

CHROMSIZES="/vast/projects/Immune_Aging/Reference_genomes/mm39/UCSC_mm39_references/mm39.chrom.sizes"

ASSEMBLY="mm39"
MIN_MAPQ=10
WALKS_POLICY="5unique"
MAX_INTER_ALIGN_GAP=20
DEDUP_MAX_MISMATCH=3
MIN_CIS_DISTANCE=1000
CPUS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p \
    "$QC_DIR/pairtools_stats" \
    "$LOG_DIR" \
    "$FULL_VALID_DIR" \
    "$TASK_TMP"

# Preserve intermediates if a command fails. Delete them only after the entire
# sample finishes successfully, preventing avoidable recomputation after an
# error late in the workflow.
SUCCESS=0
cleanup() {
    if [[ "$SUCCESS" -eq 1 ]]; then
        rm -rf "$TASK_TMP"
    else
        echo "Job failed; temporary files retained in: $TASK_TMP" >&2
    fi
}
trap cleanup EXIT

module load miniconda3
conda activate cooler-env
hash -r

for PROGRAM in samtools pairtools; do
    command -v "$PROGRAM" >/dev/null 2>&1 || {
        echo "ERROR: $PROGRAM is unavailable in cooler-env." >&2
        exit 1
    }
done

# Select one input BAM deterministically for this array task.
mapfile -t BAM_FILES < <(
    find "$BAM_DIR" -maxdepth 1 -type f \
        -name '*_Name_sorted_BL_removed.bam' -print | LC_ALL=C sort
)

N_BAMS=${#BAM_FILES[@]}
(( N_BAMS > 0 )) || { echo "ERROR: No matching BAMs in $BAM_DIR" >&2; exit 1; }
(( SLURM_ARRAY_TASK_ID < N_BAMS )) || {
    echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} exceeds $N_BAMS BAMs." >&2
    exit 1
}

BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
FILENAME=$(basename "$BAM_FILE")
SUFFIX="_Name_sorted_BL_removed.bam"
[[ "$FILENAME" == *"$SUFFIX" ]] || { echo "ERROR: Unexpected filename $FILENAME" >&2; exit 1; }
SAMPLE="${FILENAME%${SUFFIX}}"

samtools quickcheck -v "$BAM_FILE"

echo "Sample: $SAMPLE"
echo "Input BAM: $BAM_FILE"
echo "samtools: $(command -v samtools)"
echo "pairtools: $(command -v pairtools)"

PARSED_PAIRS="${TASK_TMP}/${SAMPLE}.parsed.pairs.gz"
SORTED_PAIRS="${TASK_TMP}/${SAMPLE}.sorted.pairs.gz"
NODUP_PAIRS="${TASK_TMP}/${SAMPLE}.nodup.pairs.gz"

# This full-depth valid file is needed by Script 06b to establish the common
# downsampling depth. It is deliberately not named "downsampled".
FULL_VALID_PAIRS="${FULL_VALID_DIR}/${SAMPLE}_BL_removed_deduped_valid_full.pairs.gz"

PARSE_STATS="${QC_DIR}/pairtools_stats/${SAMPLE}_parse.stats"
DEDUP_STATS="${QC_DIR}/pairtools_stats/${SAMPLE}_dedup.stats"
FULL_VALID_STATS="${QC_DIR}/pairtools_stats/${SAMPLE}_full_valid.stats"

rm -f "$PARSED_PAIRS" "$SORTED_PAIRS" "$NODUP_PAIRS" \
      "$FULL_VALID_PAIRS" "$PARSE_STATS" "$DEDUP_STATS" "$FULL_VALID_STATS"

# 1. PARSE
# Convert grouped BAM alignments into Hi-C pair records. MIN_MAPQ=10 follows
# the previous diffHic preparePairs threshold. The walk settings describe how
# pairtools resolves chimeric Hi-C alignments.
echo "[1/5] pairtools parse"
pairtools parse \
    --chroms-path "$CHROMSIZES" \
    --assembly "$ASSEMBLY" \
    --min-mapq "$MIN_MAPQ" \
    --walks-policy "$WALKS_POLICY" \
    --max-inter-align-gap "$MAX_INTER_ALIGN_GAP" \
    --nproc-in "$CPUS" \
    --nproc-out "$CPUS" \
    --output-stats "$PARSE_STATS" \
    --output "$PARSED_PAIRS" \
    "$BAM_FILE"

[[ -s "$PARSED_PAIRS" ]] || { echo "ERROR: parse produced no output." >&2; exit 1; }

# 2. SORT
# Deduplication requires upper-triangular, coordinate-sorted pairs.
echo "[2/5] pairtools sort"
pairtools sort \
    --nproc "$CPUS" \
    --tmpdir "$TASK_TMP" \
    --output "$SORTED_PAIRS" \
    "$PARSED_PAIRS"

[[ -s "$SORTED_PAIRS" ]] || { echo "ERROR: sort produced no output." >&2; exit 1; }

# 3. DEDUPLICATE
# PCR/optical duplicates and unmapped records are discarded, while their
# numbers remain available in the built-in dedup statistics file.
echo "[3/5] pairtools dedup"
pairtools dedup \
    --max-mismatch "$DEDUP_MAX_MISMATCH" \
    --output "$NODUP_PAIRS" \
    --output-dups /dev/null \
    --output-unmapped /dev/null \
    --output-stats "$DEDUP_STATS" \
    "$SORTED_PAIRS"

[[ -s "$NODUP_PAIRS" ]] || { echo "ERROR: dedup retained no mapped pairs." >&2; exit 1; }

# 4. SELECT FINAL FULL-DEPTH VALID CONTACTS
# Keep pairs that are uniquely mapped at both ends (UU). For cis pairs, remove
# contacts shorter than 1 kb; these are enriched for unligated/dangling-end and
# other very short-range artefacts. All trans contacts are retained.
echo "[4/5] pairtools select"
PAIR_SELECTION="(pair_type == 'UU') and ((chrom1 != chrom2) or (abs(pos2 - pos1) >= ${MIN_CIS_DISTANCE}))"

pairtools select \
    --nproc-in "$CPUS" \
    --nproc-out "$CPUS" \
    --output "$FULL_VALID_PAIRS" \
    "$PAIR_SELECTION" \
    "$NODUP_PAIRS"

[[ -s "$FULL_VALID_PAIRS" ]] || { echo "ERROR: select produced no valid pairs." >&2; exit 1; }

# 5. QC AT THE ACTUAL FULL-DEPTH MATRIX-INPUT STAGE
# pairtools stats supplies total, cis/trans, distance-threshold and pair-type
# metrics without loading the complete file into memory.
echo "[5/5] pairtools stats"
pairtools stats \
    --with-chromsizes \
    --n-dist-bins-decade 8 \
    --output "$FULL_VALID_STATS" \
    "$FULL_VALID_PAIRS"

FULL_VALID_COUNT=$(awk -F '\t' '$1 == "total" {print $2; exit}' "$FULL_VALID_STATS")
[[ "$FULL_VALID_COUNT" =~ ^[0-9]+$ ]] || {
    echo "ERROR: Could not read total from $FULL_VALID_STATS" >&2
    exit 1
}

echo "Full-depth valid pairs: $FULL_VALID_COUNT"
echo "Intermediate full-depth valid file: $FULL_VALID_PAIRS"
echo "Run Script 06b only after all six array tasks complete."

SUCCESS=1
