#!/bin/bash
#SBATCH --job-name=hic_downsample_cool
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/07_hic_downsample_cool_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/07_hic_downsample_cool_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

set -euo pipefail
shopt -s nullglob

# -----------------------------------------------------------------------------
# PURPOSE
# After all six samples have been parsed, deduplicated and filtered by Script
# 06, determine the smallest full-depth valid-pair library, randomly sample
# every larger library to the corresponding fraction, and construct the final
# matched-depth Cooler matrices.
#
# pairtools sample performs Bernoulli sampling by fraction, so realised counts
# will be extremely close to, but not necessarily exactly equal to, the target.
# The realised count for every sample is recorded and checked against each
# Cooler matrix sum.
# -----------------------------------------------------------------------------

PROJECT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_BL_removed_coolers_pairs_files"
OUTPUT_DIR="${PROJECT_DIR}/downsampled_merged_cool_files"
PAIRS_DIR="${PROJECT_DIR}/downsampled_merged_pairs_files"
QC_DIR="${PROJECT_DIR}/downsampled_merged_HiC_QC"
LOG_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs"

FULL_VALID_DIR="/vast/scratch/users/stacpoole.q/HiC_full_depth_valid_pairs"
TASK_TMP="/vast/scratch/users/stacpoole.q/hic_downsample_cool_${SLURM_JOB_ID}"

CHROMSIZES="/vast/projects/Immune_Aging/Reference_genomes/mm39/UCSC_mm39_references/mm39.chrom.sizes"
ASSEMBLY="mm39"
RESOLUTIONS="50000,100000,200000"
EXPECTED_SAMPLE_COUNT=6
RANDOM_SEED=20260731
CPUS="${SLURM_CPUS_PER_TASK:-8}"

# Change to "yes" only after verifying all final pairs, stats and Coolers.
# Leaving this as "no" retains full-depth valid pairs in scratch until its
# routine two-week cleanup, allowing easy reruns with different depths.
DELETE_FULL_DEPTH_VALID_AFTER_SUCCESS="no"

mkdir -p \
    "$OUTPUT_DIR" \
    "$PAIRS_DIR" \
    "$QC_DIR/pairtools_stats" \
    "$QC_DIR/combined_tables" \
    "$LOG_DIR" \
    "$TASK_TMP"

module load miniconda3
conda activate cooler-env
hash -r

for PROGRAM in pairtools cooler; do
    command -v "$PROGRAM" >/dev/null 2>&1 || {
        echo "ERROR: $PROGRAM is unavailable in cooler-env." >&2
        exit 1
    }
done

IFS=',' read -r -a RESOLUTION_ARRAY <<< "$RESOLUTIONS"
FULL_VALID_FILES=("$FULL_VALID_DIR"/*_BL_removed_deduped_valid_full.pairs.gz)

if (( ${#FULL_VALID_FILES[@]} != EXPECTED_SAMPLE_COUNT )); then
    echo "ERROR: Expected $EXPECTED_SAMPLE_COUNT full-depth valid pairs files," >&2
    echo "       but found ${#FULL_VALID_FILES[@]} in $FULL_VALID_DIR" >&2
    exit 1
fi

# Read full-depth totals from the small pairtools stats files made by Script
# 06. This avoids rescanning large pairs files merely to determine the target.
declare -A SAMPLE_FILE
declare -A SAMPLE_COUNT
TARGET_COUNT=""

COUNT_TABLE="${QC_DIR}/combined_tables/full_depth_valid_pair_counts.tsv"
printf 'sample\tfull_depth_valid_pairs\tfull_depth_pairs_file\n' > "$COUNT_TABLE"

for FULL_FILE in "${FULL_VALID_FILES[@]}"; do
    NAME=$(basename "$FULL_FILE")
    SAMPLE="${NAME%_BL_removed_deduped_valid_full.pairs.gz}"
    STATS_FILE="${QC_DIR}/pairtools_stats/${SAMPLE}_full_valid.stats"

    [[ -s "$STATS_FILE" ]] || { echo "ERROR: Missing $STATS_FILE" >&2; exit 1; }

    COUNT=$(awk -F '\t' '$1 == "total" {print $2; exit}' "$STATS_FILE")
    [[ "$COUNT" =~ ^[0-9]+$ ]] || { echo "ERROR: Invalid total in $STATS_FILE" >&2; exit 1; }

    SAMPLE_FILE["$SAMPLE"]="$FULL_FILE"
    SAMPLE_COUNT["$SAMPLE"]="$COUNT"

    if [[ -z "$TARGET_COUNT" ]] || (( COUNT < TARGET_COUNT )); then
        TARGET_COUNT="$COUNT"
    fi

    printf '%s\t%s\t%s\n' "$SAMPLE" "$COUNT" "$FULL_FILE" >> "$COUNT_TABLE"
done

echo "Smallest full-depth valid library: $TARGET_COUNT pairs"

# Record both requested and realised depths. The realised downsampled count is
# the count that must equal the matrix-wide sum in each resulting Cooler.
DOWNSAMPLE_TABLE="${QC_DIR}/combined_tables/downsampling_and_cooler_validation.tsv"
printf 'sample\tfull_depth_valid_pairs\ttarget_pairs\tsampling_fraction\trealised_downsampled_pairs\tresolution\tcooler_sum\tstatus\n' > "$DOWNSAMPLE_TABLE"

# Sort sample names for reproducible processing and table order.
mapfile -t SAMPLES < <(printf '%s\n' "${!SAMPLE_FILE[@]}" | LC_ALL=C sort)

for SAMPLE in "${SAMPLES[@]}"; do
    FULL_FILE="${SAMPLE_FILE[$SAMPLE]}"
    CURRENT_COUNT="${SAMPLE_COUNT[$SAMPLE]}"

    FINAL_PAIRS="${PAIRS_DIR}/${SAMPLE}_BL_removed_deduped_downsampled.pairs.gz"
    FINAL_STATS="${QC_DIR}/pairtools_stats/${SAMPLE}_BL_removed_deduped_downsampled.stats"

    rm -f "$FINAL_PAIRS" "$FINAL_STATS"

    if (( CURRENT_COUNT == TARGET_COUNT )); then
        # The smallest library already defines the target and is copied without
        # stochastic resampling, preserving all of its valid contacts.
        FRACTION="1.000000000000"
        cp "$FULL_FILE" "$FINAL_PAIRS"
    else
        FRACTION=$(awk -v target="$TARGET_COUNT" -v current="$CURRENT_COUNT" \
            'BEGIN {printf "%.12f", target/current}')

        pairtools sample \
            --seed "$RANDOM_SEED" \
            --nproc-in "$CPUS" \
            --nproc-out "$CPUS" \
            --output "$FINAL_PAIRS" \
            "$FRACTION" \
            "$FULL_FILE"
    fi

    [[ -s "$FINAL_PAIRS" ]] || { echo "ERROR: No final pairs for $SAMPLE" >&2; exit 1; }

    pairtools stats \
        --with-chromsizes \
        --n-dist-bins-decade 8 \
        --output "$FINAL_STATS" \
        "$FINAL_PAIRS"

    REALISED_COUNT=$(awk -F '\t' '$1 == "total" {print $2; exit}' "$FINAL_STATS")
    [[ "$REALISED_COUNT" =~ ^[0-9]+$ ]] || {
        echo "ERROR: Could not read realised total for $SAMPLE" >&2
        exit 1
    }

    echo "$SAMPLE: full=$CURRENT_COUNT target=$TARGET_COUNT realised=$REALISED_COUNT fraction=$FRACTION"

    for RAW_RESOLUTION in "${RESOLUTION_ARRAY[@]}"; do
        RESOLUTION="${RAW_RESOLUTION//[[:space:]]/}"
        COOL_FILE="${OUTPUT_DIR}/${SAMPLE}_BL_removed_deduped_downsampled_${RESOLUTION}.cool"
        COOL_INFO="${QC_DIR}/pairtools_stats/${SAMPLE}_BL_removed_deduped_downsampled_${RESOLUTION}_cooler.info"

        rm -f "$COOL_FILE" "$COOL_INFO"

        cooler cload pairs \
            --assembly "$ASSEMBLY" \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            "${CHROMSIZES}:${RESOLUTION}" \
            "$FINAL_PAIRS" \
            "$COOL_FILE"

        [[ -s "$COOL_FILE" ]] || { echo "ERROR: Missing $COOL_FILE" >&2; exit 1; }

        cooler info "$COOL_FILE" | tee "$COOL_INFO"

        # cooler info is JSON in Cooler 0.10.4. Extract the matrix-wide contact
        # sum with awk, avoiding a custom Python program.
        COOL_SUM=$(awk -F ': *' '
            /"sum"/ {
                value=$2
                gsub(/[ ,]/, "", value)
                print value
                exit
            }
        ' "$COOL_INFO")

        if [[ "$COOL_SUM" == "$REALISED_COUNT" ]]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$SAMPLE" "$CURRENT_COUNT" "$TARGET_COUNT" "$FRACTION" \
            "$REALISED_COUNT" "$RESOLUTION" "$COOL_SUM" "$STATUS" \
            >> "$DOWNSAMPLE_TABLE"

        [[ "$STATUS" == "PASS" ]] || {
            echo "ERROR: Cooler sum ($COOL_SUM) != realised pairs ($REALISED_COUNT) for $SAMPLE" >&2
            exit 1
        }
    done
done

if [[ "$DELETE_FULL_DEPTH_VALID_AFTER_SUCCESS" == "yes" ]]; then
    rm -f "${FULL_VALID_FILES[@]}"
    echo "Deleted full-depth valid-pair intermediates from scratch."
else
    echo "Full-depth valid-pair intermediates retained in $FULL_VALID_DIR"
fi

rm -rf "$TASK_TMP"

echo "Finished making matched-depth deduplicated valid pairs and Coolers."
echo "Count table: $COUNT_TABLE"
echo "Validation table: $DOWNSAMPLE_TABLE"
