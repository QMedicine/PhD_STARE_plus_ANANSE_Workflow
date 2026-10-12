#!/bin/bash

# Input directory containing TAD data
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files/hicFindTADs_TADs"

# Output directory for plots
OUTPUT_DIR="${INPUT_DIR}/OvY_TAD_comparison_metrics_plots"
mkdir -p "${OUTPUT_DIR}"

# Load miniconda and activate the hicexplorer environment
module load miniconda3
conda activate hicexplorer

# Python script for generating plots with added checks
PYTHON_PLOT_SCRIPT=$(cat <<'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Input arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]

# Font settings
#plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

# Resolutions to process
resolutions = [10000, 20000, 30000, 50000, 100000, 200000, 500000, 1000000, 2000000]

# Initialize lists to collect data across resolutions
size_data = []
insulation_data = []

# Process each resolution
for resolution in resolutions:
    # File paths for Young and Old
    old_domains_file = os.path.join(input_dir, f"O_LTHSC_{resolution}_MMP_combined_replicates_merged_hicNormalize_hicCorrectMatrix_domains.bed")
    young_domains_file = os.path.join(input_dir, f"Y_LTHSC_{resolution}_MMP_combined_replicates_merged_hicNormalize_hicCorrectMatrix_domains.bed")
    old_scores_file = os.path.join(input_dir, f"O_LTHSC_{resolution}_MMP_combined_replicates_merged_hicNormalize_hicCorrectMatrix_score.bedgraph")
    young_scores_file = os.path.join(input_dir, f"Y_LTHSC_{resolution}_MMP_combined_replicates_merged_hicNormalize_hicCorrectMatrix_score.bedgraph")
    
    # Check if files exist
    if not (os.path.isfile(old_domains_file) and os.path.isfile(young_domains_file) and 
            os.path.isfile(old_scores_file) and os.path.isfile(young_scores_file)):
        print(f"Missing files for resolution {resolution}, skipping.")
        continue
    
    print(f"\nProcessing resolution: {resolution}")
    print(f"Old domains file: {old_domains_file}")
    print(f"Young domains file: {young_domains_file}")
    print(f"Old scores file: {old_scores_file}")
    print(f"Young scores file: {young_scores_file}")
    
    # Load TAD domains and calculate sizes
    for file, condition in [(old_domains_file, "Old"), (young_domains_file, "Young")]:
        domains = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "ID", "Score", "Strand", "Extra_Start", "Extra_End", "RGB"])
        domains["Size"] = domains["End"] - domains["Start"]
        size_data.append(pd.DataFrame({"Size": domains["Size"], "Condition": condition, "Resolution": resolution}))
    
    # Load insulation scores
    for file, condition in [(old_scores_file, "Old"), (young_scores_file, "Young")]:
        scores = pd.read_csv(file, sep="\t", header=None, names=["Chromosome", "Start", "End", "Score"])
        insulation_data.append(pd.DataFrame({"Score": scores["Score"], "Condition": condition, "Resolution": resolution}))

# Combine data across resolutions
size_df = pd.concat(size_data, ignore_index=True)
insulation_df = pd.concat(insulation_data, ignore_index=True)

# Print extracted data for verification
print("\n=== Size Data Head ===")
print(size_df.head())
print("\n=== Insulation Data Head ===")
print(insulation_df.head())

# Calculate medians and quartiles for TAD sizes
print("\n=== TAD Size Statistics ===")
for resolution in resolutions:
    subset = size_df[size_df["Resolution"] == resolution]
    if not subset.empty:
        print(f"Resolution: {resolution}")
        print(f"  Median: {subset['Size'].median()}")
        print(f"  Q1: {subset['Size'].quantile(0.25)}")
        print(f"  Q3: {subset['Size'].quantile(0.75)}")

# Calculate medians and quartiles for insulation scores
print("\n=== Insulation Score Statistics ===")
for resolution in resolutions:
    subset = insulation_df[insulation_df["Resolution"] == resolution]
    if not subset.empty:
        print(f"Resolution: {resolution}")
        print(f"  Median: {subset['Score'].median()}")
        print(f"  Q1: {subset['Score'].quantile(0.25)}")
        print(f"  Q3: {subset['Score'].quantile(0.75)}")

# Create violin plot for TAD sizes
sns.violinplot(data=size_df, x="Resolution", y="Size", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"}, split=True)
plt.title("TAD Sizes by Resolution", fontsize=14)
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("TAD Size", fontsize=12)
plt.legend(title="Condition", fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Size_violin_plot.pdf"))
plt.close()

# Create violin plot for insulation scores
sns.violinplot(data=insulation_df, x="Resolution", y="Score", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"}, split=True)
plt.title("TAD Insulation Scores by Resolution", fontsize=14)
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("Insulation Score", fontsize=12)
plt.legend(title="Condition", fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Insulation_violin_plot.pdf"))
plt.close()
EOF
)

# Run the Python script
python3 -c "${PYTHON_PLOT_SCRIPT}" "${INPUT_DIR}" "${OUTPUT_DIR}"

# Deactivate the conda environment
conda deactivate

echo "Plots saved to ${OUTPUT_DIR}"
