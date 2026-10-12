#!/bin/bash

# Input directory containing TAD results
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files/hicFindTADs_TADs"

# Output directory for results
OUTPUT_DIR="${INPUT_DIR}/OvY_TAD_comparison_metrics"
mkdir -p "${OUTPUT_DIR}"

# Python script for analysis with added checks
PYTHON_SCRIPT=$(cat <<'EOF'
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import sys

# Input arguments
old_domains_file = sys.argv[1]
young_domains_file = sys.argv[2]
old_scores_file = sys.argv[3]
young_scores_file = sys.argv[4]
output_file = sys.argv[5]

# Load domain files and ensure numeric types
def load_domains(file):
    domains = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "ID", "Score", "Strand", "Extra_Start", "Extra_End", "RGB"])
    domains["Start"] = pd.to_numeric(domains["Start"], errors="coerce")
    domains["End"] = pd.to_numeric(domains["End"], errors="coerce")
    domains.dropna(inplace=True)  # Remove rows with non-numeric values
    domains["Size"] = domains["End"] - domains["Start"]
    return domains[["Chromosome", "Start", "End", "Size"]]  # Keep only relevant columns

# Load score files and ensure numeric types
def load_scores(file):
    scores = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "Score"])
    scores["Score"] = pd.to_numeric(scores["Score"], errors="coerce")
    scores.dropna(inplace=True)  # Remove rows with non-numeric values
    return scores

# Load data
old_domains = load_domains(old_domains_file)
young_domains = load_domains(young_domains_file)
old_scores = load_scores(old_scores_file)
young_scores = load_scores(young_scores_file)

# Output data summaries for debugging
print("=== Old Domains Data ===")
print(old_domains.head())
print(old_domains.describe())
print("\n")

print("=== Young Domains Data ===")
print(young_domains.head())
print(young_domains.describe())
print("\n")

print("=== Old Scores Data ===")
print(old_scores.head())
print(old_scores.describe())
print("\n")

print("=== Young Scores Data ===")
print(young_scores.head())
print(young_scores.describe())
print("\n")

# Check for empty datasets
if old_domains.empty or young_domains.empty or old_scores.empty or young_scores.empty:
    print(f"Empty dataset detected. Skipping resolution.")
    with open(output_file, "w") as f:
        f.write("Empty dataset detected. Skipping resolution.\n")
    sys.exit(0)

# Perform calculations
results = {
    "Metric": ["Total TADs", "Median Size", "Mean Size", "Median Insulation", "Mean Insulation", "Size Significance", "Insulation Significance"],
    "Old": [
        len(old_domains),
        np.median(old_domains["Size"]),
        np.mean(old_domains["Size"]),
        np.median(old_scores["Score"]),
        np.mean(old_scores["Score"]),
        None,  # Placeholder for significance
        None   # Placeholder for significance
    ],
    "Young": [
        len(young_domains),
        np.median(young_domains["Size"]),
        np.mean(young_domains["Size"]),
        np.median(young_scores["Score"]),
        np.mean(young_scores["Score"]),
        None,  # Placeholder for significance
        None   # Placeholder for significance
    ]
}

# Perform significance tests
_, results["Old"][5] = mannwhitneyu(old_domains["Size"], young_domains["Size"])
_, results["Old"][6] = mannwhitneyu(old_scores["Score"], young_scores["Score"])

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)
EOF
)

# Loop through resolutions
for resolution in 10000 20000 30000 50000 100000 200000 500000 1000000 2000000; do
    echo "Processing resolution: ${resolution}"

    # Input files for old and young conditions
    OLD_DOMAINS="${INPUT_DIR}/O_LTHSC_${resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_domains.bed"
    YOUNG_DOMAINS="${INPUT_DIR}/Y_LTHSC_${resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_domains.bed"
    OLD_SCORES="${INPUT_DIR}/O_LTHSC_${resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_score.bedgraph"
    YOUNG_SCORES="${INPUT_DIR}/Y_LTHSC_${resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_score.bedgraph"

    # Output file for results
    OUTPUT_FILE="${OUTPUT_DIR}/OvY_TAD_comparison_metrics_${resolution}.csv"

    # Check if files exist
    if [[ -f "${OLD_DOMAINS}" && -f "${YOUNG_DOMAINS}" && -f "${OLD_SCORES}" && -f "${YOUNG_SCORES}" ]]; then
        # Run Python analysis
        python3 -c "${PYTHON_SCRIPT}" "${OLD_DOMAINS}" "${YOUNG_DOMAINS}" "${OLD_SCORES}" "${YOUNG_SCORES}" "${OUTPUT_FILE}"
        echo "Results saved to: ${OUTPUT_FILE}"
    else
        echo "Missing files for resolution: ${resolution}. Skipping."
    fi

    echo "--------------------------------------------"
done

echo "TAD comparison analysis complete."
