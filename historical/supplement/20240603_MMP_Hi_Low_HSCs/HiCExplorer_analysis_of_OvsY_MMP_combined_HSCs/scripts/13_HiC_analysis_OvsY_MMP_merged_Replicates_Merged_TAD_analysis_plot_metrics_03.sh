#!/bin/bash

# Input directory containing TAD data
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files/hicFindTADs_TADs"

# Output directory for plots
OUTPUT_DIR="${INPUT_DIR}/OvY_TAD_comparison_metrics/OvY_TAD_comparison_metrics_plots"
mkdir -p "${OUTPUT_DIR}"

# Load miniconda and activate the hicexplorer environment
module load miniconda3
conda activate hicexplorer

# Python script for generating box plots
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
plt.rcParams['font.size'] = 12

# Resolutions to process
resolutions = [10000, 20000, 30000, 50000, 100000, 200000]

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

# Create box plot for TAD sizes
sns.boxplot(data=size_df, x="Resolution", y="Size", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"})
plt.title("TAD Sizes by Resolution", fontsize=14)
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("TAD Size", fontsize=12)
plt.ylim(0, 0.75 * 1e7)
plt.legend(title="Condition", fontsize=12, loc='upper right', bbox_to_anchor=(1.30, 1))
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Size_box_plot.pdf"))
plt.close()

# Create box plot for insulation scores
sns.boxplot(data=insulation_df, x="Resolution", y="Score", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"})
plt.title("TAD Insulation Scores by Resolution", fontsize=14)
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("Insulation Score", fontsize=12)
plt.ylim(-5.0, 5.0)
plt.legend(title="Condition", fontsize=12, loc='upper right', bbox_to_anchor=(1.30, 1))
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Insulation_box_plot.pdf"))
plt.close()
EOF
)

# Run the Python script
python3 -c "${PYTHON_PLOT_SCRIPT}" "${INPUT_DIR}" "${OUTPUT_DIR}"

# Deactivate the conda environment
conda deactivate

echo "Box plots saved to ${OUTPUT_DIR}"
