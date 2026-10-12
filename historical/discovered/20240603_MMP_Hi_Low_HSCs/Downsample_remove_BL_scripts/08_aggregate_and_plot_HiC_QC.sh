#!/bin/bash
#SBATCH --job-name=aggregate_hic_qc
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/08_aggregate_hic_qc_%j.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs/08_aggregate_hic_qc_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
shopt -s nullglob

# =============================================================================
# PURPOSE
# =============================================================================
# Aggregate the small QC files produced by:
#   06_convert_bams_to_pairs.sh
#   07_downsample_valid_pairs_and_create_coolers.sh
#
# This script does not read BAM files or load complete pairs files into memory.
# It uses pairtools statistics, the downsampling validation table and metadata
# stored in the final Cooler files. The resulting tables and vector PDF plots
# compare all six libraries side-by-side.
# =============================================================================

# =============================================================================
# Paths
# =============================================================================

PROJECT_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/downsampled_merged_BL_removed_coolers_pairs_files"

OUTPUT_DIR="${PROJECT_DIR}/downsampled_merged_cool_files"
PAIRS_DIR="${PROJECT_DIR}/downsampled_merged_pairs_files"
QC_DIR="${PROJECT_DIR}/downsampled_merged_HiC_QC"

STATS_DIR="${QC_DIR}/pairtools_stats"
SOURCE_TABLE_DIR="${QC_DIR}/combined_tables"
TABLE_DIR="${QC_DIR}/combined_tables/08_aggregate_QC"
PLOT_DIR="${QC_DIR}/combined_plots/08_aggregate_QC"

LOG_DIR="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Logs/bam_to_cool_logs"

EXPECTED_SAMPLE_COUNT=6

# Leave empty to use alphabetic sample order. To specify the order used in all
# plots, provide comma-separated sample names exactly as they appear in files:
# SAMPLE_ORDER="Young_1,Young_2,Young_3,Old_1,Old_2,Old_3"
SAMPLE_ORDER=""

mkdir -p \
    "$TABLE_DIR" \
    "$PLOT_DIR" \
    "$LOG_DIR"

# =============================================================================
# Software
# =============================================================================

module load miniconda3
conda activate cooler-env
hash -r

for PROGRAM in python pairtools cooler; do
    command -v "$PROGRAM" >/dev/null 2>&1 || {
        echo "ERROR: $PROGRAM is unavailable in cooler-env." >&2
        exit 1
    }
done

echo "============================================================"
echo "Aggregating matched-depth Hi-C QC"
echo "Date:       $(date)"
echo "Python:     $(command -v python)"
echo "Pairtools:  $(command -v pairtools)"
echo "Cooler:     $(command -v cooler)"
echo "============================================================"

# =============================================================================
# Lightweight input checks
# =============================================================================

PARSE_STATS=("$STATS_DIR"/*_parse.stats)
DEDUP_STATS=("$STATS_DIR"/*_dedup.stats)
FULL_VALID_STATS=("$STATS_DIR"/*_full_valid.stats)
DOWNSAMPLED_STATS=("$STATS_DIR"/*_BL_removed_deduped_downsampled.stats)
FINAL_PAIRS=("$PAIRS_DIR"/*_BL_removed_deduped_downsampled.pairs.gz)
COOL_FILES=("$OUTPUT_DIR"/*_BL_removed_deduped_downsampled_*.cool)

COUNT_TABLE="${SOURCE_TABLE_DIR}/full_depth_valid_pair_counts.tsv"
DOWNSAMPLE_TABLE="${SOURCE_TABLE_DIR}/downsampling_and_cooler_validation.tsv"

for DESCRIPTION in \
    "parse stats:${#PARSE_STATS[@]}" \
    "dedup stats:${#DEDUP_STATS[@]}" \
    "full-depth valid stats:${#FULL_VALID_STATS[@]}" \
    "downsampled stats:${#DOWNSAMPLED_STATS[@]}" \
    "final downsampled pairs:${#FINAL_PAIRS[@]}"
do
    NAME=${DESCRIPTION%%:*}
    COUNT=${DESCRIPTION##*:}

    if (( COUNT != EXPECTED_SAMPLE_COUNT )); then
        echo "ERROR: Expected ${EXPECTED_SAMPLE_COUNT} ${NAME} files but found ${COUNT}." >&2
        exit 1
    fi
done

if (( ${#COOL_FILES[@]} == 0 )); then
    echo "ERROR: No final Cooler files were found in $OUTPUT_DIR" >&2
    exit 1
fi

[[ -s "$COUNT_TABLE" ]] || {
    echo "ERROR: Missing full-depth count table: $COUNT_TABLE" >&2
    exit 1
}

[[ -s "$DOWNSAMPLE_TABLE" ]] || {
    echo "ERROR: Missing downsampling validation table: $DOWNSAMPLE_TABLE" >&2
    exit 1
}

# =============================================================================
# Aggregate and plot
# =============================================================================
# Only small text tables and Cooler metadata are read. The final pairs files are
# checked for existence but are not rescanned, because Script 07 already checked
# their pairtools totals against the Cooler sums.
# =============================================================================

python - \
    "$STATS_DIR" \
    "$PAIRS_DIR" \
    "$OUTPUT_DIR" \
    "$COUNT_TABLE" \
    "$DOWNSAMPLE_TABLE" \
    "$TABLE_DIR" \
    "$PLOT_DIR" \
    "$EXPECTED_SAMPLE_COUNT" \
    "$SAMPLE_ORDER" <<'PY'
import math
import re
import sys
from pathlib import Path

import cooler
import matplotlib

matplotlib.use("pdf")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

(
    stats_dir,
    pairs_dir,
    cooler_dir,
    count_table_path,
    downsample_table_path,
    table_dir,
    plot_dir,
    expected_sample_count,
    requested_sample_order,
) = sys.argv[1:]

stats_dir = Path(stats_dir)
pairs_dir = Path(pairs_dir)
cooler_dir = Path(cooler_dir)
count_table_path = Path(count_table_path)
downsample_table_path = Path(downsample_table_path)
table_dir = Path(table_dir)
plot_dir = Path(plot_dir)
expected_sample_count = int(expected_sample_count)

table_dir.mkdir(parents=True, exist_ok=True)
plot_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Publication-friendly vector PDF settings.
# Type 42 embedding keeps text selectable/editable in Affinity Designer.
# -----------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Arial",
            "Helvetica",
            "Liberation Sans",
            "DejaVu Sans",
        ],
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 9.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    }
)


def clean_axes(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def read_pairtools_stats(path):
    """Read the key/value rows in a pairtools stats file."""
    values = {}

    with open(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                continue

            key = fields[0]
            text = fields[1]

            try:
                if text.lower() in {"nan", ".nan"}:
                    value = math.nan
                elif re.fullmatch(r"[-+]?\d+", text):
                    value = int(text)
                else:
                    value = float(text)
            except ValueError:
                value = text

            values[key] = value

    return values


def require_integer(stats, key, path):
    if key not in stats:
        raise RuntimeError(f"Required key '{key}' is absent from {path}")

    value = stats[key]
    if isinstance(value, float) and math.isnan(value):
        raise RuntimeError(f"Key '{key}' is NaN in {path}")

    return int(value)


def stat_integer(stats, key, default=0):
    value = stats.get(key, default)
    if isinstance(value, float) and math.isnan(value):
        return int(default)
    return int(value)


def derive_sample(filename, suffix):
    if not filename.endswith(suffix):
        raise RuntimeError(f"Unexpected filename: {filename}")
    return filename[: -len(suffix)]


# -----------------------------------------------------------------------------
# Determine samples from the final downsampled stats files.
# -----------------------------------------------------------------------------
final_stats_suffix = "_BL_removed_deduped_downsampled.stats"
final_stats_paths = sorted(stats_dir.glob(f"*{final_stats_suffix}"))

samples_found = [
    derive_sample(path.name, final_stats_suffix)
    for path in final_stats_paths
]

if len(samples_found) != expected_sample_count:
    raise RuntimeError(
        f"Expected {expected_sample_count} samples but found "
        f"{len(samples_found)}: {samples_found}"
    )

if len(samples_found) != len(set(samples_found)):
    raise RuntimeError(f"Duplicate sample names were found: {samples_found}")

if requested_sample_order.strip():
    samples = [x.strip() for x in requested_sample_order.split(",") if x.strip()]
    if set(samples) != set(samples_found) or len(samples) != len(samples_found):
        raise RuntimeError(
            "SAMPLE_ORDER must contain every sample exactly once.\n"
            f"Found samples: {samples_found}\n"
            f"Requested: {samples}"
        )
else:
    samples = sorted(samples_found)

print("Sample order:")
for sample in samples:
    print(f"  {sample}")

# -----------------------------------------------------------------------------
# Read the two compact tables written by Script 07.
# -----------------------------------------------------------------------------
full_depth_table = pd.read_csv(count_table_path, sep="\t")
downsample_table = pd.read_csv(downsample_table_path, sep="\t")

if set(full_depth_table["sample"]) != set(samples):
    raise RuntimeError("Samples in full_depth_valid_pair_counts.tsv do not match stats files.")

if set(downsample_table["sample"]) != set(samples):
    raise RuntimeError("Samples in downsampling_and_cooler_validation.tsv do not match stats files.")

# One downsampling row per sample is sufficient for sample-level quantities;
# Script 07 repeats these values once per Cooler resolution.
downsample_per_sample = (
    downsample_table.sort_values(["sample", "resolution"])
    .drop_duplicates("sample")
    .set_index("sample")
)

# -----------------------------------------------------------------------------
# Build one QC record per sample from pairtools' built-in statistics.
# -----------------------------------------------------------------------------
summary_rows = []

for sample in samples:
    parse_path = stats_dir / f"{sample}_parse.stats"
    dedup_path = stats_dir / f"{sample}_dedup.stats"
    full_valid_path = stats_dir / f"{sample}_full_valid.stats"
    final_path = stats_dir / f"{sample}{final_stats_suffix}"
    final_pairs_path = pairs_dir / f"{sample}_BL_removed_deduped_downsampled.pairs.gz"

    for required_path in [parse_path, dedup_path, full_valid_path, final_path, final_pairs_path]:
        if not required_path.exists():
            raise RuntimeError(f"Missing required file: {required_path}")

    parse_stats = read_pairtools_stats(parse_path)
    dedup_stats = read_pairtools_stats(dedup_path)
    full_stats = read_pairtools_stats(full_valid_path)
    final_stats = read_pairtools_stats(final_path)

    # The dedup statistics provide the clearest mutually related processing
    # counts. Some pairtools versions omit a few optional keys, so defaults are
    # used only for non-essential subcategories.
    total_pairs = require_integer(dedup_stats, "total", dedup_path)
    unmapped = stat_integer(dedup_stats, "total_unmapped")
    single_sided = stat_integer(dedup_stats, "total_single_sided_mapped")
    duplicates = stat_integer(dedup_stats, "total_dups")
    nodups = require_integer(dedup_stats, "total_nodups", dedup_path)

    full_valid = require_integer(full_stats, "total", full_valid_path)
    downsampled_valid = require_integer(final_stats, "total", final_path)

    # Contacts removed by pairtools select after deduplication are mostly
    # non-UU pairs and cis contacts shorter than 1 kb.
    removed_by_select = nodups - full_valid
    if removed_by_select < 0:
        raise RuntimeError(
            f"{sample}: full valid count ({full_valid}) exceeds total_nodups ({nodups})."
        )

    # Account for any pair category not explicitly represented by the common
    # pairtools keys. Usually this is zero; retaining it prevents misleading
    # stacked bars if a version reports additional categories.
    represented = unmapped + single_sided + duplicates + removed_by_select + full_valid
    other_pairs = total_pairs - represented
    if other_pairs < 0:
        # Pairtools totals can differ in scope between versions. Do not silently
        # draw a negative category; record the discrepancy for inspection.
        other_pairs = 0
        composition_difference = represented - total_pairs
    else:
        composition_difference = 0

    full_cis = stat_integer(full_stats, "cis")
    full_trans = stat_integer(full_stats, "trans")
    final_cis = stat_integer(final_stats, "cis")
    final_trans = stat_integer(final_stats, "trans")

    if final_cis + final_trans != downsampled_valid:
        raise RuntimeError(
            f"{sample}: final cis + trans does not equal final total "
            f"({final_cis} + {final_trans} != {downsampled_valid})."
        )

    final_cis_1kb = stat_integer(final_stats, "cis_1kb+", final_cis)
    final_cis_10kb = stat_integer(final_stats, "cis_10kb+")
    final_cis_1_to_10kb = max(final_cis_1kb - final_cis_10kb, 0)

    target_count = int(downsample_per_sample.loc[sample, "target_pairs"])
    table_realised = int(downsample_per_sample.loc[sample, "realised_downsampled_pairs"])
    sampling_fraction = float(downsample_per_sample.loc[sample, "sampling_fraction"])

    if table_realised != downsampled_valid:
        raise RuntimeError(
            f"{sample}: Script 07 table reports {table_realised} realised pairs, "
            f"but final pairtools stats report {downsampled_valid}."
        )

    mapped_for_duplication = duplicates + nodups
    duplicate_percentage_mapped = (
        100.0 * duplicates / mapped_for_duplication
        if mapped_for_duplication > 0
        else 0.0
    )

    summary_rows.append(
        {
            "sample": sample,
            "total_parsed_pairs": total_pairs,
            "unmapped_pairs": unmapped,
            "single_sided_mapped_pairs": single_sided,
            "duplicate_pairs": duplicates,
            "deduplicated_pairs_before_select": nodups,
            "removed_by_final_select": removed_by_select,
            "full_depth_valid_pairs": full_valid,
            "target_pairs": target_count,
            "sampling_fraction": sampling_fraction,
            "downsampled_valid_pairs": downsampled_valid,
            "downsampling_difference_from_target": downsampled_valid - target_count,
            "downsampling_difference_percent": (
                100.0 * (downsampled_valid - target_count) / target_count
                if target_count > 0
                else 0.0
            ),
            "other_or_unaccounted_pairs": other_pairs,
            "composition_count_discrepancy": composition_difference,
            "duplicate_percentage_of_mapped_pairs": duplicate_percentage_mapped,
            "full_valid_percentage_of_parsed_pairs": (
                100.0 * full_valid / total_pairs if total_pairs > 0 else 0.0
            ),
            "downsampled_percentage_of_full_valid": (
                100.0 * downsampled_valid / full_valid if full_valid > 0 else 0.0
            ),
            "full_valid_cis": full_cis,
            "full_valid_trans": full_trans,
            "downsampled_valid_cis": final_cis,
            "downsampled_valid_trans": final_trans,
            "downsampled_cis_1_to_10kb": final_cis_1_to_10kb,
            "downsampled_cis_over_10kb": final_cis_10kb,
            "downsampled_cis_percentage": (
                100.0 * final_cis / downsampled_valid if downsampled_valid > 0 else 0.0
            ),
            "downsampled_trans_percentage": (
                100.0 * final_trans / downsampled_valid if downsampled_valid > 0 else 0.0
            ),
        }
    )

summary = pd.DataFrame(summary_rows)
summary["sample"] = pd.Categorical(summary["sample"], categories=samples, ordered=True)
summary = summary.sort_values("sample").reset_index(drop=True)
summary["sample"] = summary["sample"].astype(str)

summary.to_csv(table_dir / "all_samples_HiC_QC_summary.tsv", sep="\t", index=False)

# -----------------------------------------------------------------------------
# Read metadata from the actual final Coolers and independently compare their
# matrix-wide sums with the final downsampled pairtools totals.
# -----------------------------------------------------------------------------
cooler_rows = []

for sample in samples:
    expected_count = int(
        summary.loc[summary["sample"].eq(sample), "downsampled_valid_pairs"].iloc[0]
    )

    sample_coolers = sorted(
        cooler_dir.glob(f"{sample}_BL_removed_deduped_downsampled_*.cool")
    )
    if not sample_coolers:
        raise RuntimeError(f"No final Cooler files found for {sample}")

    for cool_path in sample_coolers:
        match = re.search(r"_BL_removed_deduped_downsampled_(\d+)\.cool$", cool_path.name)
        if match is None:
            raise RuntimeError(f"Cannot determine resolution from {cool_path.name}")

        filename_resolution = int(match.group(1))
        clr = cooler.Cooler(str(cool_path))
        stored_resolution = int(clr.binsize)
        cooler_sum = int(clr.info["sum"])
        nnz = int(clr.info["nnz"])
        nbins = int(clr.info["nbins"])

        # Cooler stores the upper triangle. This occupancy measure is useful
        # for comparing matrix sparsity at the same resolution.
        possible_upper_pixels = nbins * (nbins + 1) / 2
        upper_triangle_occupancy = (
            100.0 * nnz / possible_upper_pixels if possible_upper_pixels > 0 else 0.0
        )

        cooler_rows.append(
            {
                "sample": sample,
                "resolution": stored_resolution,
                "filename_resolution": filename_resolution,
                "cool_file": str(cool_path),
                "expected_downsampled_pairs": expected_count,
                "cooler_sum": cooler_sum,
                "nnz": nnz,
                "nbins": nbins,
                "upper_triangle_occupancy_percent": upper_triangle_occupancy,
                "count_match": cooler_sum == expected_count,
                "resolution_match": stored_resolution == filename_resolution,
            }
        )

cooler_summary = pd.DataFrame(cooler_rows)
cooler_summary["sample"] = pd.Categorical(
    cooler_summary["sample"], categories=samples, ordered=True
)
cooler_summary = cooler_summary.sort_values(["sample", "resolution"]).reset_index(drop=True)
cooler_summary["sample"] = cooler_summary["sample"].astype(str)

cooler_summary.to_csv(
    table_dir / "all_samples_final_cooler_validation.tsv", sep="\t", index=False
)

if not cooler_summary["count_match"].all():
    failures = cooler_summary.loc[~cooler_summary["count_match"]]
    raise RuntimeError(
        "Some final Cooler sums differ from the final pairtools totals:\n"
        + failures.to_string(index=False)
    )

if not cooler_summary["resolution_match"].all():
    failures = cooler_summary.loc[~cooler_summary["resolution_match"]]
    raise RuntimeError(
        "Some stored Cooler resolutions differ from their filenames:\n"
        + failures.to_string(index=False)
    )

# Save a concise Methods/QC table with the most important quantities.
methods_table = summary[
    [
        "sample",
        "total_parsed_pairs",
        "duplicate_pairs",
        "duplicate_percentage_of_mapped_pairs",
        "full_depth_valid_pairs",
        "full_valid_percentage_of_parsed_pairs",
        "target_pairs",
        "sampling_fraction",
        "downsampled_valid_pairs",
        "downsampling_difference_from_target",
        "downsampling_difference_percent",
        "downsampled_valid_cis",
        "downsampled_valid_trans",
        "downsampled_cis_percentage",
        "downsampled_trans_percentage",
    ]
]
methods_table.to_csv(
    table_dir / "all_samples_methods_ready_QC_summary.tsv", sep="\t", index=False
)

# -----------------------------------------------------------------------------
# Plot helpers and common coordinates.
# -----------------------------------------------------------------------------
x = np.arange(len(samples))

# Plot 1: full processing composition before downsampling.
# This shows where pairs were lost during mapping, deduplication and selection.
composition_columns = [
    "unmapped_pairs",
    "single_sided_mapped_pairs",
    "duplicate_pairs",
    "removed_by_final_select",
    "full_depth_valid_pairs",
    "other_or_unaccounted_pairs",
]
composition_labels = [
    "Unmapped",
    "Single-sided mapped",
    "Duplicate",
    "Removed by final select",
    "Full-depth valid",
    "Other/unaccounted",
]

# Hide the final residual category when it is zero for every sample.
if summary["other_or_unaccounted_pairs"].sum() == 0:
    composition_columns = composition_columns[:-1]
    composition_labels = composition_labels[:-1]

fig, axis = plt.subplots(figsize=(6.6, 3.8))
bottom = np.zeros(len(samples))
for column, label in zip(composition_columns, composition_labels):
    values = summary[column].to_numpy()
    axis.bar(x, values, bottom=bottom, width=0.72, label=label)
    bottom += values

axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Read pairs")
axis.set_title("Hi-C processing yield: absolute counts")
axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_processing_composition_counts.pdf")
plt.close(fig)

# Plot 2: the same composition as percentages of parsed pairs.
fig, axis = plt.subplots(figsize=(6.6, 3.8))
bottom = np.zeros(len(samples))
for column, label in zip(composition_columns, composition_labels):
    values = 100.0 * summary[column].to_numpy() / summary["total_parsed_pairs"].to_numpy()
    axis.bar(x, values, bottom=bottom, width=0.72, label=label)
    bottom += values

axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Parsed read pairs (%)")
axis.set_ylim(0, max(100, float(np.nanmax(bottom)) * 1.02))
axis.set_title("Hi-C processing yield: proportions")
axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_processing_composition_proportions.pdf")
plt.close(fig)

# Plot 3: duplication and usable full-depth valid-pair yield.
width = 0.36
fig, axis = plt.subplots(figsize=(5.3, 3.3))
axis.bar(
    x - width / 2,
    summary["duplicate_percentage_of_mapped_pairs"],
    width=width,
    label="Duplicates among mapped pairs",
)
axis.bar(
    x + width / 2,
    summary["full_valid_percentage_of_parsed_pairs"],
    width=width,
    label="Full-depth valid pairs",
)
axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Read pairs (%)")
axis.set_title("Duplicate fraction and final usable yield")
axis.legend(frameon=False)
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_duplicate_and_valid_yield.pdf")
plt.close(fig)

# Plot 4: full-depth valid counts, target, and realised downsampled counts.
fig, axis = plt.subplots(figsize=(6.0, 3.5))
axis.bar(x, summary["full_depth_valid_pairs"], width=0.64, label="Full-depth valid")
axis.scatter(x, summary["target_pairs"], marker="_", s=180, linewidths=2.0, label="Target")
axis.scatter(x, summary["downsampled_valid_pairs"], marker="x", s=38, linewidths=1.2, label="Realised")
axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Valid contacts")
axis.set_title("Matched-depth downsampling")
axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axis.legend(frameon=False)
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_downsampling_depths.pdf")
plt.close(fig)

# Plot 5: stochastic deviation from the target after Bernoulli sampling.
fig, axis = plt.subplots(figsize=(5.1, 3.1))
axis.axhline(0, linewidth=0.7)
axis.bar(x, summary["downsampling_difference_percent"], width=0.64)
axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Difference from target (%)")
axis.set_title("Realised downsampling depth relative to target")
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_downsampling_target_deviation.pdf")
plt.close(fig)

# Plot 6: cis/trans proportions in the final matrix-input pairs.
final_cis_pct = summary["downsampled_cis_percentage"].to_numpy()
final_trans_pct = summary["downsampled_trans_percentage"].to_numpy()
fig, axis = plt.subplots(figsize=(5.2, 3.3))
axis.bar(x, final_cis_pct, width=0.7, label="Cis")
axis.bar(x, final_trans_pct, bottom=final_cis_pct, width=0.7, label="Trans")
axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Downsampled valid pairs (%)")
axis.set_ylim(0, 100)
axis.set_title("Final cis and trans contact proportions")
axis.legend(frameon=False)
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_final_cis_trans.pdf")
plt.close(fig)

# Plot 7: final valid-pair genomic spans available from pairtools stats.
span_columns = [
    "downsampled_cis_1_to_10kb",
    "downsampled_cis_over_10kb",
    "downsampled_valid_trans",
]
span_labels = ["Cis 1–10 kb", "Cis >10 kb", "Trans"]
fig, axis = plt.subplots(figsize=(5.5, 3.5))
bottom = np.zeros(len(samples))
for column, label in zip(span_columns, span_labels):
    values = 100.0 * summary[column].to_numpy() / summary["downsampled_valid_pairs"].to_numpy()
    axis.bar(x, values, bottom=bottom, width=0.7, label=label)
    bottom += values

axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Downsampled valid pairs (%)")
axis.set_ylim(0, 100)
axis.set_title("Final valid-pair genomic spans")
axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_final_pair_spans.pdf")
plt.close(fig)

# Plot 8: number of occupied Cooler pixels by matrix resolution.
# At matched contact depth, this provides a simple view of relative sparsity.
fig, axis = plt.subplots(figsize=(5.8, 3.5))
for resolution in sorted(cooler_summary["resolution"].unique()):
    subset = (
        cooler_summary.loc[cooler_summary["resolution"].eq(resolution)]
        .set_index("sample")
        .reindex(samples)
    )
    axis.plot(x, subset["nnz"], marker="o", label=f"{resolution:,} bp")

axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Non-zero matrix pixels")
axis.set_title("Final Cooler matrix occupancy")
axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axis.legend(frameon=False)
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_cooler_nonzero_pixels.pdf")
plt.close(fig)

# Plot 9: final pairtools totals and all Cooler sums. Symbols should overlap.
fig, axis = plt.subplots(figsize=(6.0, 3.5))
resolutions = sorted(cooler_summary["resolution"].unique())
offsets = np.linspace(-0.22, 0.22, len(resolutions))
for offset, resolution in zip(offsets, resolutions):
    subset = (
        cooler_summary.loc[cooler_summary["resolution"].eq(resolution)]
        .set_index("sample")
        .reindex(samples)
    )
    axis.scatter(
        x + offset,
        subset["cooler_sum"],
        s=22,
        label=f"{resolution:,} bp Cooler",
    )

axis.scatter(
    x,
    summary["downsampled_valid_pairs"],
    marker="x",
    s=45,
    linewidths=1.2,
    label="Final pairtools total",
)
axis.set_xticks(x)
axis.set_xticklabels(samples, rotation=30, ha="right")
axis.set_ylabel("Contact count")
axis.set_title("Final pairs-to-Cooler count validation")
axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axis.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
clean_axes(axis)
fig.savefig(plot_dir / "all_samples_pairs_to_cooler_validation.pdf")
plt.close(fig)

print("\nAll final Cooler sums and resolutions passed validation.\n")
print(summary.to_string(index=False))
print("\nCooler metadata:\n")
print(cooler_summary.to_string(index=False))
print(f"\nTables: {table_dir}")
print(f"Plots:  {plot_dir}")
PY

echo
echo "============================================================"
echo "Combined Hi-C QC completed successfully"
echo "Tables: $TABLE_DIR"
echo "Plots:  $PLOT_DIR"
echo "Completed: $(date)"
echo "============================================================"