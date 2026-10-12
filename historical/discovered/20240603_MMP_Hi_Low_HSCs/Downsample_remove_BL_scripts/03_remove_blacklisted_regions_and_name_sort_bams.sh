#!/bin/bash
#SBATCH --job-name=remove_blacklist
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/03_remove_blacklisted_regions_and_name_sort_bams/03_remove_blacklisted_regions_and_name_sort_bams_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/03_remove_blacklisted_regions_and_name_sort_bams/03_remove_blacklisted_regions_and_name_sort_bams_%A_%a.err
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

set -euo pipefail
shopt -s nullglob

# ---------------------------------------------------------------------------
# Directories and files
# ---------------------------------------------------------------------------

BAM_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files"

BLACKLIST="/vast/projects/Immune_Aging/Reference_genomes/mm39/mm39.excluderanges.bed"

OUTPUT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/COSorted_merged_bam_files_BL_removed"

OUTPUT_DIR2="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Name_sorted_merged_bam_files_BL_removed"

LOG_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/03_remove_blacklisted_regions_and_name_sort_bams"

SCRATCH_ROOT="/vast/scratch/users/stacpoole.q"

mkdir -p \
    "$OUTPUT_DIR" \
    "$OUTPUT_DIR2" \
    "$LOG_DIR" \
    "$SCRATCH_ROOT"

# Give every Slurm array task its own temporary directory.
TASK_TMP="${SCRATCH_ROOT}/remove_blacklist_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$TASK_TMP"

# Delete this task's temporary directory when the job finishes or fails.
trap 'rm -rf "$TASK_TMP"' EXIT

# ---------------------------------------------------------------------------
# Software
# ---------------------------------------------------------------------------

module load samtools
module load bedtools

CPUS="${SLURM_CPUS_PER_TASK:-4}"

echo "============================================================"
echo "Blacklist filtering started"
echo "Date:                 $(date)"
echo "Job ID:               ${SLURM_JOB_ID}"
echo "Array task ID:        ${SLURM_ARRAY_TASK_ID}"
echo "Allocated CPUs:       ${CPUS}"
echo "Temporary directory:  ${TASK_TMP}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Identify the BAM assigned to this array task
# ---------------------------------------------------------------------------

mapfile -t BAM_FILES < <(
    find "$BAM_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*_COsorted.bam' \
        -print |
    LC_ALL=C sort
)

N_BAMS=${#BAM_FILES[@]}

if (( N_BAMS == 0 )); then
    echo "ERROR: No files matching *_COsorted.bam were found in:"
    echo "       $BAM_DIR"
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID >= N_BAMS )); then
    echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} is outside the BAM list."
    echo "       Number of input BAM files: ${N_BAMS}"
    exit 1
fi

BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
FILENAME=$(basename "$BAM_FILE")

EXPECTED_SUFFIX="_COsorted.bam"

if [[ "$FILENAME" != *"$EXPECTED_SUFFIX" ]]; then
    echo "ERROR: Unexpected input filename:"
    echo "       $FILENAME"
    echo "Expected suffix:"
    echo "       $EXPECTED_SUFFIX"
    exit 1
fi

SAMPLE="${FILENAME%${EXPECTED_SUFFIX}}"

# ---------------------------------------------------------------------------
# Output and temporary filenames
# ---------------------------------------------------------------------------

FILTERED_BAM="${OUTPUT_DIR}/${SAMPLE}_COsorted_BL_removed.bam"

FINAL_BAM="${OUTPUT_DIR2}/${SAMPLE}_Name_sorted_BL_removed.bam"

BLACKLIST_QNAMES="${TASK_TMP}/${SAMPLE}.blacklist_qnames.txt"

OVERLAP_CHECK_BAM="${TASK_TMP}/${SAMPLE}.retained_blacklist_overlap.bam"

echo
echo "Input BAM:"
echo "  $BAM_FILE"
echo
echo "Blacklist BED:"
echo "  $BLACKLIST"
echo
echo "Coordinate-sorted output:"
echo "  $FILTERED_BAM"
echo
echo "Coordinate-sorted BAM index:"
echo "  ${FILTERED_BAM}.bai"
echo
echo "Name-sorted output:"
echo "  $FINAL_BAM"
echo

# ---------------------------------------------------------------------------
# Preliminary validation
# ---------------------------------------------------------------------------

if [[ ! -s "$BAM_FILE" ]]; then
    echo "ERROR: Input BAM does not exist or is empty:"
    echo "       $BAM_FILE"
    exit 1
fi

if [[ ! -s "$BLACKLIST" ]]; then
    echo "ERROR: Blacklist BED does not exist or is empty:"
    echo "       $BLACKLIST"
    exit 1
fi

if ! samtools quickcheck -v "$BAM_FILE"; then
    echo "ERROR: samtools quickcheck failed for:"
    echo "       $BAM_FILE"
    exit 1
fi

# ---------------------------------------------------------------------------
# Count input records
# ---------------------------------------------------------------------------

echo "Counting records in the input BAM..."

BEFORE_RECORDS=$(samtools view \
    -@ "$CPUS" \
    -c \
    "$BAM_FILE")

echo "Input alignment records: $BEFORE_RECORDS"

# ---------------------------------------------------------------------------
# Identify all QNAMEs for which at least one alignment overlaps a blacklist
# interval.
# ---------------------------------------------------------------------------

echo
echo "Identifying read-pair names associated with blacklist overlaps..."

bedtools intersect \
    -u \
    -abam "$BAM_FILE" \
    -b "$BLACKLIST" |
samtools view - |
cut -f 1 |
LC_ALL=C sort \
    --unique \
    --temporary-directory="$TASK_TMP" \
    > "$BLACKLIST_QNAMES"

BLACKLIST_QNAME_COUNT=$(wc -l < "$BLACKLIST_QNAMES")
BLACKLIST_QNAME_COUNT=${BLACKLIST_QNAME_COUNT//[[:space:]]/}

echo "Unique query names marked for pair-level removal: $BLACKLIST_QNAME_COUNT"

# ---------------------------------------------------------------------------
# Remove every alignment associated with a blacklist-overlapping QNAME.
#
# Because the input BAM is coordinate sorted and filtering does not reorder
# records, the resulting FILTERED_BAM remains coordinate sorted.
# ---------------------------------------------------------------------------

echo
echo "Removing complete read pairs..."

if (( BLACKLIST_QNAME_COUNT > 0 )); then
    samtools view \
        -@ "$CPUS" \
        -b \
        -N "^${BLACKLIST_QNAMES}" \
        -o "$FILTERED_BAM" \
        "$BAM_FILE"
else
    echo "No blacklist-overlapping query names were found."
    echo "Writing the input alignments unchanged to the output BAM."

    samtools view \
        -@ "$CPUS" \
        -b \
        -o "$FILTERED_BAM" \
        "$BAM_FILE"
fi

if ! samtools quickcheck -v "$FILTERED_BAM"; then
    echo "ERROR: samtools quickcheck failed for the filtered BAM:"
    echo "       $FILTERED_BAM"
    exit 1
fi

# ---------------------------------------------------------------------------
# Index the coordinate-sorted filtered BAM
# ---------------------------------------------------------------------------

echo
echo "Indexing the coordinate-sorted filtered BAM..."

# Remove an old index first so a failed indexing attempt cannot leave a stale
# index appearing to belong to the new BAM.
rm -f \
    "${FILTERED_BAM}.bai" \
    "${FILTERED_BAM%.bam}.bai" \
    "${FILTERED_BAM}.csi"

samtools index \
    -@ "$CPUS" \
    "$FILTERED_BAM"

if [[ -s "${FILTERED_BAM}.bai" ]]; then
    FILTERED_INDEX="${FILTERED_BAM}.bai"
elif [[ -s "${FILTERED_BAM%.bam}.bai" ]]; then
    FILTERED_INDEX="${FILTERED_BAM%.bam}.bai"
elif [[ -s "${FILTERED_BAM}.csi" ]]; then
    FILTERED_INDEX="${FILTERED_BAM}.csi"
else
    echo "ERROR: No BAM index was produced for:"
    echo "       $FILTERED_BAM"
    exit 1
fi

echo "PASS: Coordinate-sorted BAM index created:"
echo "      $FILTERED_INDEX"

# Confirm that samtools can use the index to obtain index statistics.
if ! samtools idxstats "$FILTERED_BAM" > "${TASK_TMP}/${SAMPLE}.idxstats.txt"; then
    echo "ERROR: samtools idxstats could not read the filtered BAM index."
    exit 1
fi

echo "PASS: samtools idxstats successfully read the BAM index."

# ---------------------------------------------------------------------------
# Produce the requested name-sorted BAM
# ---------------------------------------------------------------------------

echo
echo "Name-sorting the filtered BAM..."

samtools sort \
    -n \
    -@ "$CPUS" \
    -T "${TASK_TMP}/${SAMPLE}.namesort" \
    -o "$FINAL_BAM" \
    "$FILTERED_BAM"

if ! samtools quickcheck -v "$FINAL_BAM"; then
    echo "ERROR: samtools quickcheck failed for the name-sorted BAM:"
    echo "       $FINAL_BAM"
    exit 1
fi

# Do not run samtools index on FINAL_BAM. It is name sorted, not coordinate
# sorted, and therefore is not suitable for a standard BAM genomic index.

# ---------------------------------------------------------------------------
# Validation 1: compare record counts before and after filtering
# ---------------------------------------------------------------------------

echo
echo "Counting records after pair-level blacklist removal..."

AFTER_RECORDS=$(samtools view \
    -@ "$CPUS" \
    -c \
    "$FILTERED_BAM")

NAME_SORTED_RECORDS=$(samtools view \
    -@ "$CPUS" \
    -c \
    "$FINAL_BAM")

REMOVED_RECORDS=$(( BEFORE_RECORDS - AFTER_RECORDS ))

if (( BEFORE_RECORDS > 0 )); then
    RETAINED_PERCENT=$(awk \
        -v after="$AFTER_RECORDS" \
        -v before="$BEFORE_RECORDS" \
        'BEGIN {printf "%.4f", (after / before) * 100}')

    REMOVED_PERCENT=$(awk \
        -v removed="$REMOVED_RECORDS" \
        -v before="$BEFORE_RECORDS" \
        'BEGIN {printf "%.4f", (removed / before) * 100}')
else
    RETAINED_PERCENT="0.0000"
    REMOVED_PERCENT="0.0000"
fi

echo
echo "Record-count validation"
echo "-----------------------"
echo "Input records:                    $BEFORE_RECORDS"
echo "Coordinate-sorted output records: $AFTER_RECORDS"
echo "Name-sorted output records:       $NAME_SORTED_RECORDS"
echo "Removed records:                  $REMOVED_RECORDS"
echo "Retained percentage:              ${RETAINED_PERCENT}%"
echo "Removed percentage:               ${REMOVED_PERCENT}%"

if (( AFTER_RECORDS != NAME_SORTED_RECORDS )); then
    echo
    echo "ERROR: Record count changed during name sorting."
    echo "Coordinate-sorted BAM: $AFTER_RECORDS"
    echo "Name-sorted BAM:       $NAME_SORTED_RECORDS"
    exit 1
else
    echo "PASS: Coordinate-sorted and name-sorted BAM record counts match."
fi

# ---------------------------------------------------------------------------
# Validation 2: verify that no retained alignment overlaps the blacklist
# ---------------------------------------------------------------------------

echo
echo "Checking retained alignments for blacklist overlap..."

bedtools intersect \
    -u \
    -abam "$FILTERED_BAM" \
    -b "$BLACKLIST" \
    > "$OVERLAP_CHECK_BAM"

RETAINED_OVERLAP_RECORDS=$(samtools view \
    -c \
    "$OVERLAP_CHECK_BAM")

echo "Retained records overlapping blacklist: $RETAINED_OVERLAP_RECORDS"

if (( RETAINED_OVERLAP_RECORDS != 0 )); then
    echo "FAIL: The filtered BAM still contains alignments overlapping the blacklist."
    exit 1
else
    echo "PASS: No retained alignment overlaps a blacklist interval."
fi

# ---------------------------------------------------------------------------
# Validation 3: verify that no excluded QNAME remains
# ---------------------------------------------------------------------------

echo
echo "Checking that no records from excluded query names remain..."

if (( BLACKLIST_QNAME_COUNT > 0 )); then
    RETAINED_EXCLUDED_QNAME_RECORDS=$(samtools view \
        -@ "$CPUS" \
        -c \
        -N "$BLACKLIST_QNAMES" \
        "$FILTERED_BAM")
else
    RETAINED_EXCLUDED_QNAME_RECORDS=0
fi

echo "Retained records carrying an excluded query name: $RETAINED_EXCLUDED_QNAME_RECORDS"

if (( RETAINED_EXCLUDED_QNAME_RECORDS != 0 )); then
    echo "FAIL: At least one record belonging to a blacklisted read pair remains."
    exit 1
else
    echo "PASS: All records belonging to blacklist-associated query names were removed."
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

echo
echo "============================================================"
echo "Blacklist filtering completed successfully"
echo "Sample: $SAMPLE"
echo
echo "Pair-level filtering:"
echo "  Unique query names excluded: $BLACKLIST_QNAME_COUNT"
echo "  Alignment records before:    $BEFORE_RECORDS"
echo "  Alignment records after:     $AFTER_RECORDS"
echo "  Alignment records removed:   $REMOVED_RECORDS"
echo "  Retained percentage:         ${RETAINED_PERCENT}%"
echo "  Removed percentage:          ${REMOVED_PERCENT}%"
echo
echo "Validation:"
echo "  Coordinate/name-sort record counts match: PASS"
echo "  Retained blacklist overlaps:              0"
echo "  Retained excluded-QNAME records:          0"
echo "  Coordinate-sorted BAM index usable:       PASS"
echo
echo "Outputs:"
echo "  Coordinate-sorted BAM:"
echo "    $FILTERED_BAM"
echo "  Coordinate-sorted BAM index:"
echo "    $FILTERED_INDEX"
echo "  Name-sorted BAM:"
echo "    $FINAL_BAM"
echo
echo "Temporary files were stored in:"
echo "  $TASK_TMP"
echo
echo "Temporary files will now be removed."
echo "Date completed: $(date)"
echo "============================================================"