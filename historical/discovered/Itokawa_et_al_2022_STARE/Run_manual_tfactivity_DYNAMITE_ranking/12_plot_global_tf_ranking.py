#!/usr/bin/env python3
# 12_plot_global_tf_ranking.py
import argparse, pandas as pd, numpy as np
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

p = argparse.ArgumentParser()
p.add_argument("--tf_rank_annot", required=True)   # ...tf_ranking.annot.tsv from step 9
p.add_argument("--coefs", required=True)           # ...ATAC-seq.filtered.tsv  (TF,value)
p.add_argument("--topN", type=int, default=20)
p.add_argument("--flip-coef", action="store_true")
p.add_argument("--summary", help="...topTF_topGenes.summary.with_gabc.tsv")
p.add_argument("--out", required=True)
args = p.parse_args()

rank = pd.read_csv(args.tf_rank_annot, sep="\t", index_col=0)  # index = tf_col; has 'dcg','jaspar_id','tf_symbol_clean'
coefs = pd.read_csv(args.coefs, sep="\t", index_col=0)         # index = tf_col; col 'value'
if args.flip_coef:
    coefs["value"] = -coefs["value"]

df = rank.copy()
df["coef"] = coefs.reindex(df.index)["value"].astype(float)
df["negcoef"] = df["coef"] < 0

tf_de = None
if args.summary:
    S = pd.read_csv(args.summary, sep="\t")
    if {"tf","TF_log2FC","TF_padj"}.issubset(S.columns):
        tf_de = (S.drop_duplicates(subset=["tf"])
                   .set_index("tf")[["TF_log2FC","TF_padj"]])
        df = df.join(tf_de, how="left")

def color_from_de(lfc, padj, alpha=0.05):
    if not pd.notnull(padj) or padj >= alpha:
        return (0.65,0.65,0.65)  # grey
    return (0.85,0.2,0.2) if (lfc or 0) > 0 else (0.2,0.35,0.85)

top = df.sort_values("dcg", ascending=False).head(args.topN).copy()
colors = [color_from_de(top.loc[i,"TF_log2FC"], top.loc[i,"TF_padj"]) for i in top.index]

# Build display labels: "TF [ID]" (use the exact motif id; if multiple semicolon IDs were ever present, take first)
def pick_label(row):
    sym = str(row.get("tf_symbol_clean") or row.name)
    jid = str(row.get("jaspar_id") or "").split(";")[0]
    return f"{sym} [{jid}]" if jid else sym

labels = top.apply(pick_label, axis=1)

fig, ax = plt.subplots(figsize=(10,5))
bars = ax.bar(labels, top["dcg"], color=colors, edgecolor="black", linewidth=0.4)
for bar, is_neg in zip(bars, top["negcoef"].tolist()):
    if is_neg:
        bar.set_hatch("//")

ax.set_ylabel("DCG (1 − (rank−1)/N)")
ax.set_title("Top TFs by DCG\n(hatch = negative coefficient;\ncolor = TF DE: red↑old, blue↓old, grey ns)")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")  # << fix for 'ha' error
fig.tight_layout()
fig.savefig(args.out, dpi=300)
