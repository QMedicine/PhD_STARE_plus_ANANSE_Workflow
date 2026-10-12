#!/bin/bash

# Input directory containing TAD data
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files/hicCorrectMatrix_corrected_h5_files/hicFindTADs_TADs"

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
resolutions = [50000]

# Initialize lists to collect data across resolutions
size_data = []
insulation_data = []

# Process each resolution
for resolution in resolutions:
    # File paths for Young and Old
    old_domains_file = os.path.join(input_dir, f"O_LTHSC_{resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_domains.bed")
    young_domains_file = os.path.join(input_dir, f"Y_LTHSC_{resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_domains.bed")
    old_scores_file = os.path.join(input_dir, f"O_LTHSC_{resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_score.bedgraph")
    young_scores_file = os.path.join(input_dir, f"Y_LTHSC_{resolution}_MMP_combined_replicates_merged_BL_removed_hicNormalize_hicCorrectMatrix_score.bedgraph")
    
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

#############################
# P-values for the comparisons
size_pval = 0.060986406915173874
insulation_pval = 5.120467764747628e-06

# Function to determine significance stars
def get_significance_stars(pval):
    if pval <= 0.0001:
        return '****'
    elif pval <= 0.001:
        return '***'
    elif pval <= 0.01:
        return '**'
    elif pval <= 0.05:
        return '*'
    else:
        return 'ns'

# Get stars for the p-values
size_stars = get_significance_stars(size_pval)
insulation_stars = get_significance_stars(insulation_pval)

# Positions for significance stars and connecting lines
size_y_position = 0.28 * 1e7  # Adjusted y-position for TAD Size
size_line_height = size_y_position - 0.05 * 1e7  # Connecting line height

insulation_y_position = 1.8  # Adjusted y-position for Insulation Score
insulation_line_height = insulation_y_position - 0.2  # Connecting line height

##############################

# Set the width of the figure
fig_width = 3  # Narrower figure width
fig_height = 4

# Create box plot for TAD sizes (without dots)
plt.figure(figsize=(fig_width, fig_height))
sns.boxplot(data=size_df, x="Resolution", y="Size", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"}, showfliers=False)
plt.title("Old vs Young HSC TAD\nSizes (50 kbp resolution)", fontsize=14, loc='center')
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("TAD Size", fontsize=12)
plt.ylim(0, size_y_position * 1.1)  # Slightly increase ylim to accommodate stars
plt.xticks(rotation=45)  # Rotate x-axis labels by 45 degrees
plt.legend(title="Condition", fontsize=12, loc='upper right', bbox_to_anchor=(1.40, 1))
plt.text(0.5, size_y_position, size_stars, ha='center', fontsize=14, color='black')  # Add stars

# Add significance stars and connecting line
x_coords = [0, 1]  # Indices of the two box plots
plt.plot(x_coords, [size_line_height] * 2, color='black', linewidth=1.2)  # Horizontal line
plt.plot([x_coords[0], x_coords[0]], [size_line_height, size_y_position], color='black', linewidth=1.2)  # Left vertical line
plt.plot([x_coords[1], x_coords[1]], [size_line_height, size_y_position], color='black', linewidth=1.2)  # Right vertical line
plt.text(0.5, size_y_position, size_stars, ha='center', fontsize=14, color='black')  # Add stars

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Size_box_plot_no_dots_50kbp.pdf"))
plt.close()

# Create box plot for insulation scores (without dots)
plt.figure(figsize=(fig_width, fig_height))
sns.boxplot(data=insulation_df, x="Resolution", y="Score", hue="Condition", palette={"Young": "lightgreen", "Old": "orange"}, showfliers=False)
plt.title("Old vs Young HSC TAD Insulation\nScores (50 kbp resolution)", fontsize=14, loc='center')
plt.xlabel("Resolution (bp)", fontsize=12)
plt.ylabel("Insulation Score", fontsize=12)
plt.ylim(-2.0, insulation_y_position * 1.1)  # Slightly increase ylim to accommodate stars
plt.xticks(rotation=45)  # Rotate x-axis labels by 45 degrees
plt.legend(title="Condition", fontsize=12, loc='upper right', bbox_to_anchor=(1.40, 1))
plt.text(0.5, insulation_y_position, insulation_stars, ha='center', fontsize=14, color='black')  # Add stars

# Add significance stars and connecting line
x_coords = [0, 1]  # Indices of the two box plots
plt.plot(x_coords, [insulation_line_height] * 2, color='black', linewidth=1.2)  # Horizontal line
plt.plot([x_coords[0], x_coords[0]], [insulation_line_height, insulation_y_position], color='black', linewidth=1.2)  # Left vertical line
plt.plot([x_coords[1], x_coords[1]], [insulation_line_height, insulation_y_position], color='black', linewidth=1.2)  # Right vertical line
plt.text(0.5, insulation_y_position, insulation_stars, ha='center', fontsize=14, color='black')  # Add stars


plt.tight_layout()
plt.savefig(os.path.join(output_dir, "TAD_Insulation_box_plot_no_dots_50kbp.pdf"))
plt.close()
EOF
)

# Run the Python script
python3 -c "${PYTHON_PLOT_SCRIPT}" "${INPUT_DIR}" "${OUTPUT_DIR}"

# Deactivate the conda environment
conda deactivate

echo "Box plots (without dots) saved to ${OUTPUT_DIR}"
