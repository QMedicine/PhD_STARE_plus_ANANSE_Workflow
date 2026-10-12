#!/usr/bin/env python3
# 13_plot_tf_gene_heatmap.py
import argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--summary", required=True)   # ...topTF_topGenes.summary.with_gabc.tsv (step 11 output)
p.add_argument("--tf_rank_annot", required=True)
p.add_argument("--topT", type=int, default=10)
p.add_argument("--topG", type=int, default=20)
p.add_argument("--out", required=True)       # pdf
args = p.parse_args()

summ = pd.read_csv(args.summary, sep="\t")
rank = pd.read_csv(args.tf_rank_annot, sep="\t", index_col=0)

# choose TFs
tf_order = rank.sort_values("dcg", ascending=False).head(args.topT).index.tolist()
summ = summ[summ["tf"].isin(tf_order)]

# pick topG genes per TF by score
top_genes = []
for tf in tf_order:
    sg = summ[summ["tf"]==tf].nlargest(args.topG, "score")["gene"].tolist()
    top_genes.extend(sg)
top_genes = list(dict.fromkeys(top_genes))  # preserve order, dedupe

# build matrix
M = pd.pivot_table(summ[summ["gene"].isin(top_genes)], index="gene", columns="tf", values="score", aggfunc="max").fillna(0.0)
# keep order
M = M.loc[[g for g in top_genes if g in M.index], [t for t in tf_order if t in M.columns]]

plt.figure(figsize=(0.2*len(M.columns)+2, 0.2*len(M.index)+2))
im = plt.imshow(M.values, aspect="auto", interpolation="nearest")
plt.colorbar(im, label="TF–TG score")
plt.xticks(range(M.shape[1]), M.columns, rotation=45, ha="right")
plt.yticks(range(M.shape[0]), M.index)
plt.title("Top genes (rows) per top TF (cols) by TF–TG score")
plt.tight_layout()
plt.savefig(args.out)
