# 08_create_ranking.py
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, statistics as st
from scipy import stats
# --- imports these for extra dataframe and plots ---
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

p = argparse.ArgumentParser()
p.add_argument("--scores", required=True)  # TF–TG score matrix (rows=genes, cols=TFs)
p.add_argument("--tfout", required=True)   # TF ranking (DCG per TF)
p.add_argument("--tgout", required=True)   # Gene-wise DCG per TF
a = p.parse_args()

# Get output directory for saving new files
out_dir = os.path.dirname(a.tfout)
# Use the --tfout name as a base for the new files
base_name = os.path.basename(a.tfout).split('.')[0] 

# Define output file paths
full_stats_path = os.path.join(out_dir, f"{base_name}_full_stats.tsv")
plot_path = os.path.join(out_dir, f"{base_name}_top10_boxplot.pdf")

df = pd.read_csv(a.scores, sep="\t", index_col=0)

def mwu(bg, fg):
    # Added error handling for constant arrays (a common issue)
    try:
        _, p = stats.mannwhitneyu(bg, fg, alternative='two-sided')
    except ValueError:
        p = 1.0
    return p

rank = pd.DataFrame(index=df.columns, columns=["sum","mean","q95","q99","median","bkg_median","pvalue"])
rank["sum"]   = df.sum()
rank["mean"]  = df.mean()
rank["q95"]   = df.quantile(0.95)
rank["q99"]   = df.quantile(0.99)
rank["median"]= df.median()

# --- Dynamic "Leave-One-Column-Out" Background ---
print("Calculating dynamic backgrounds and p-values...")
pvals = []
bkg_meds = []
all_cols = set(df.columns)

for col_name in df.columns:
    # Create background by dropping the current column
    bg_cols = list(all_cols - {col_name})
    bg_mat = df[bg_cols]
    
    # Flatten to a vector
    bg_vec = bg_mat.values.ravel()
    
    # Calculate and store stats
    bkg_meds.append(st.median(bg_vec))
    pvals.append(mwu(bg_vec, df[col_name].values))

rank["bkg_median"] = bkg_meds
rank["pvalue"] = pvals
# --- End of dynamic background block ---

# Calculate BH-adjusted p-values
reject, p_adjusted, _, _ = multipletests(rank["pvalue"], alpha=0.05, method='fdr_bh')
rank["p_adjusted_BH"] = p_adjusted

# Save the full, unfiltered stats table as *full_stats.tsv file which should have a different bg_median for each TF
print(f"Saving full stats table to {full_stats_path}...")
rank.sort_values(by="p_adjusted_BH", ascending=True).to_csv(full_stats_path, sep="\t")

# Filter using *dynamic* background median and *adjusted* p-value
rank = rank[(rank["median"] > rank["bkg_median"]) & (rank["p_adjusted_BH"] < 0.05)]
rank = rank.sort_values(by="median", ascending=False)

# --- Plotting Top 10 TFs vs Background ---

print(f"Generating plot and saving to {plot_path}...")

# Get the top 10 TFs from the *filtered and sorted* list
top_10_tfs = rank.head(10)
top_10_tf_names = top_10_tfs.index.tolist()

# Helper function for p-value asterisks
def p_to_asterisks(p):
    if p > 0.05: return 'ns'
    if p <= 0.0001: return '****'
    if p <= 0.001: return '***'
    if p <= 0.01: return '**'
    if p <= 0.05: return '*'
    return ''

# --- Prepare data for side-by-side plotting ---
plot_data = []
plot_order = []
color_palette = {}
all_cols = set(df.columns)

for tf_name in top_10_tf_names:
    bg_name = f"{tf_name}_BG"
    
    # Define plot order and colors
    plot_order.extend([bg_name, tf_name])
    color_palette[bg_name] = "yellow"
    color_palette[tf_name] = "pink"

    # Get TF foreground data
    for score in df[tf_name]:
        plot_data.append({"group": tf_name, "score": score})
        
    # Get TF *specific* background data
    bg_cols = list(all_cols - {tf_name})
    bg_vec = df[bg_cols].values.ravel() # Re-calculate background vector
    for score in bg_vec:
        plot_data.append({"group": bg_name, "score": score})

plot_df = pd.DataFrame(plot_data)

# --- End of new data prep ---
    
# Create the plot (adjust size as required)
plt.figure(figsize=(8, 3)) 

ax = sns.boxplot(
    data=plot_df, 
    x="group", 
    y="score", 
    order=plot_order, 
    palette=color_palette,
    showfliers=False  # Hide outliers to "un-squash" the plot
)

# Get the top of the y-axis (now set to the whiskers)
y_whisker_max = ax.get_ylim()[1] 
# Set the y-position for the text slightly above the whiskers
plot_text_y = y_whisker_max * 1.02
# Set the new y-axis limit 10% higher to make room for the text
plot_frame_top = y_whisker_max * 1.10
ax.set_ylim(None, plot_frame_top) 

# Get p-values for the top 10 TFs
p_values = top_10_tfs["p_adjusted_BH"]

# Add asterisks (now positioned over the TF box, at odd-numbered positions)
for i, tf_name in enumerate(top_10_tf_names):
    # x-position for TF is 1, 3, 5, ...
    x_pos = (i * 2) + 1 
    p_val = p_values.get(tf_name, 1.0)
    asterisks = p_to_asterisks(p_val)
    
    ax.text(x_pos, plot_text_y, asterisks, 
            horizontalalignment='center', size='medium', color='black')

# Formatting
ax.set_title("Top 10 Significant TFs vs. Their Specific Background")
ax.set_ylabel("TF-TG Score")
ax.set_xlabel("TF / Group")
# Rotate x-axis labels (font size smaller to fit)
plt.xticks(rotation=90, ha='center', fontsize=8) 
plt.tight_layout() 

# Save as text-searchable PDF
plt.savefig(plot_path, format="pdf", transparent=True)
plt.close() # Close the plot to free memory

# --- End of plotting snippet ---

# DCG per TF (same construction as tfactivity)
N = len(rank)
rank["rank"] = range(1, N+1)
rank["dcg"]  = 1 - (rank["rank"] - 1) / max(N,1)
rank[["dcg"]].to_csv(a.tfout, sep="\t")

# Gene-wise DCG per TF (only significant TFs)
sig_tfs = rank.index
if len(sig_tfs) > 0:
    df2 = df[sig_tfs]
    # For each TF, rank genes by descending score; convert to DCG = 1 - rank/num_genes
    gene_dcg = 1 - (df2.rank(ascending=False).astype(int) / len(df2.index))
    gene_dcg.to_csv(a.tgout, sep="\t")
else:
    # Write empty but valid files
    pd.DataFrame(columns=["dcg"]).to_csv(a.tfout, sep="\t")
    pd.DataFrame().to_csv(a.tgout, sep="\t")
