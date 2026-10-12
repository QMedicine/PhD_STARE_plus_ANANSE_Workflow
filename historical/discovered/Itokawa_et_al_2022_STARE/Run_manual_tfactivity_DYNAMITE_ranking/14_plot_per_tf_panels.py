#!/usr/bin/env python3
# 14_plot_per_tf_panels.py
import argparse, numpy as np, pandas as pd, matplotlib.pyplot as plt, re
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

p = argparse.ArgumentParser()
p.add_argument("--summary", required=True)
p.add_argument("--tf_rank_annot", required=True)
p.add_argument("--topTF", type=int, default=6)
p.add_argument("--topG", type=int, default=25)
p.add_argument("--flip-lfc", action="store_true")
p.add_argument("--flip-coef", action="store_true")
p.add_argument("--label-top", type=int, default=10)
p.add_argument("--label-metric", choices=["score","abs_lfc","abs_aff","abs_gabc"], default="score")
p.add_argument("--tf_tpm", required=True)
p.add_argument("--deseq", required=True)
p.add_argument("--idmap", default=None)
p.add_argument("--outdir", required=True)
args = p.parse_args()

# Data
S = pd.read_csv(args.summary, sep="\t")
R = pd.read_csv(args.tf_rank_annot, sep="\t", index_col=0)
order = R.sort_values("dcg", ascending=False).head(args.topTF).index.tolist()
S = S[S["tf"].isin(order)].copy()

# Orient signs as requested
S["log2FC_oriented"] = -S["log2FC"] if args.flip_lfc else S["log2FC"]
S["coef_oriented"]   = -S["coef_value"] if args.flip_coef else S["coef_value"]

# Label selection helper
def pick_labels(df, k, metric):
    if metric == "score":
        idx = df["score"].nlargest(k).index
    elif metric == "abs_lfc":
        idx = df["log2FC_oriented"].abs().nlargest(k).index
    elif metric == "abs_aff":
        idx = df["log2_aff_ratio"].abs().nlargest(k).index
    else:  # abs_gabc
        idx = df["log2_gabc_ratio"].abs().nlargest(k).index
    return set(idx)

# TPM table
tpm_df = pd.read_csv(args.tf_tpm, sep="\t")
tpm_df["__gene_key__"] = tpm_df.iloc[:,0].astype(str).str.upper()
tpm_df = tpm_df.set_index("__gene_key__", drop=False)
expr_cols = [c for c in tpm_df.columns if c not in ("__gene_key__", tpm_df.columns[0])]

def pick_cols(cols, pat):
    rx = re.compile(pat, re.IGNORECASE)
    return [c for c in cols if rx.search(c)]
old_cols = pick_cols(expr_cols, r"old|aged")
young_cols = pick_cols(expr_cols, r"young|yng")

# DE table (for TF DE)
deseq_df = pd.read_csv(args.deseq, sep="\t", index_col=0)
deseq_df.index = deseq_df.index.astype(str).str.split(".").str[0]
sym2ens = {}
if args.idmap:
    idmap = pd.read_csv(args.idmap, sep="\t", index_col=0)
    idmap.index = idmap.index.astype(str).str.split(".").str[0]
    idmap["gene_key"] = idmap["gene_name"].astype(str).str.upper()
    sym2ens = idmap.reset_index().set_index("gene_key").iloc[:,0].to_dict()

def tf_de_lookup(tf_symbol_upper):
    # try symbol as Ensembl or symbol column
    if tf_symbol_upper in deseq_df.index:
        row = deseq_df.loc[tf_symbol_upper]
        return float(row.get("log2FoldChange", np.nan)), float(row.get("padj", np.nan))
    ens = sym2ens.get(tf_symbol_upper)
    if ens and ens in deseq_df.index:
        row = deseq_df.loc[ens]
        return float(row.get("log2FoldChange", np.nan)), float(row.get("padj", np.nan))
    if "gene_name" in deseq_df.columns:
        hits = deseq_df.index[deseq_df["gene_name"].astype(str).str.upper() == tf_symbol_upper]
        if len(hits)>0:
            row = deseq_df.loc[hits[0]]
            return float(row.get("log2FoldChange", np.nan)), float(row.get("padj", np.nan))
    return np.nan, np.nan

def tf_tpm_groups(tf_symbol_upper):
    if tf_symbol_upper not in tpm_df.index:
        return [], []
    row = tpm_df.loc[tf_symbol_upper]
    return list(pd.to_numeric(row[old_cols], errors="coerce").dropna().values), \
           list(pd.to_numeric(row[young_cols], errors="coerce").dropna().values)

# Diverging cmap (blue–grey–red)
cmap_bgr = LinearSegmentedColormap.from_list(
    "blue_grey_red", [(0.0, "#2b6cb0"), (0.5, "#bfbfbf"), (1.0, "#c53030")]
)

import matplotlib.gridspec as gridspec
import os
os.makedirs(args.outdir, exist_ok=True)

for tf in order:
    sub = S[S["tf"]==tf].nlargest(args.topG, "score").copy()
    if sub.empty:
        continue

    tf_sym = (sub.iloc[0]['tf_symbol'] or tf)
    tf_key = str(tf_sym).upper()

    # TF-level DE/TPM
    old_vals, young_vals = tf_tpm_groups(tf_key)
    tf_lfc, tf_padj = tf_de_lookup(tf_key)

    # Set up figure
    fig = plt.figure(figsize=(13.5, 5.2))
    gs = gridspec.GridSpec(1, 4, width_ratios=[0.7, 0.7, 3.3, 4.3], wspace=0.5)

    # (0) Thin TF TPM boxplot
    ax0a = fig.add_subplot(gs[0,0])
    if old_vals or young_vals:
        ax0a.boxplot([young_vals, old_vals], labels=["Young","Old"], vert=True, widths=0.6)
    ax0a.set_title("TF TPM", fontsize=12)
    ax0a.tick_params(axis='x', labelrotation=45, labelsize=10)
    ax0a.tick_params(axis='y', labelsize=10)

    # (1) TF log2FC bar, centered & colored by DE significance
    ax0b = fig.add_subplot(gs[0,1])
    ax0b.axhline(0, color="k", linewidth=1)
    col = (0.65,0.65,0.65)
    alpha = 0.05
    if pd.notnull(tf_padj) and tf_padj < alpha:
        col = (0.85,0.2,0.2) if (tf_lfc or 0) > 0 else (0.2,0.35,0.85)
    ax0b.bar([0], [tf_lfc if pd.notnull(tf_lfc) else 0], width=0.6,
             color=col, edgecolor="black", linewidth=0.4)
    label = f"padj={tf_padj:.2e}" if pd.notnull(tf_padj) else "padj=NA"
    ytxt = (tf_lfc if pd.notnull(tf_lfc) else 0)
    ytxt += 0.08 * (1 if (tf_lfc or 0) >= 0 else -1)
    ax0b.text(0, ytxt, label, ha="center",
              va="bottom" if (tf_lfc or 0) >= 0 else "top", fontsize=10)
    ax0b.set_xticks([0]); ax0b.set_xticklabels(["TF log2FC"], rotation=45, fontsize=10)
    yl = max(1.0, abs(tf_lfc) if pd.notnull(tf_lfc) else 1.0)
    ax0b.set_ylim(-(yl * 1.1), (yl * 1.1))
    ax0b.tick_params(axis='y', labelsize=10)
    ax0b.set_title("TF DE", fontsize=12)

    # (2) Top genes (barh), colored by oriented log2FC with 0-centered diverging norm
    ax1 = fig.add_subplot(gs[0,2])
    gs_sub = sub.sort_values("score", ascending=True).reset_index(drop=True)
    vals = gs_sub["log2FC_oriented"].astype(float)
    vmin, vmax = np.nanmin(vals), np.nanmax(vals)
    if not np.isfinite(vmin): vmin = -1.0
    if not np.isfinite(vmax): vmax =  1.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    colors = cmap_bgr(norm(vals.values))
    ax1.barh(range(len(gs_sub)), gs_sub["score"], color=colors, edgecolor="none")
    ax1.set_yticks(range(len(gs_sub))); ax1.set_yticklabels(gs_sub["gene"])
    ax1.set_xlabel("TF–TG score")
    ax1.set_title(f"{tf_sym} ({tf}) — Top {len(gs_sub)} genes\n(color = oriented log2FC)")
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap_bgr),
                        ax=ax1, fraction=0.046, pad=0.01)
    cbar.set_label("oriented log2FC (Old up positive)", fontsize=10)

    # (3) Scatter (size ~ score), labels on chosen points
    ax2 = fig.add_subplot(gs[0,3])
    sizes = np.clip(gs_sub["score"]/np.nanmax(gs_sub["score"])*200, 20, 220)
    sc = ax2.scatter(gs_sub["log2_aff_ratio"], gs_sub["log2FC_oriented"],
                     s=sizes, c=gs_sub["log2FC_oriented"], cmap=cmap_bgr, norm=norm,
                     edgecolors="k", linewidths=0.25)
    ax2.axvline(0, linewidth=1); ax2.axhline(0, linewidth=1)
    ax2.set_xlabel("log2 affinity ratio (Old/Young)")
    ax2.set_ylabel("oriented log2FC (Old up positive)")
    ax2.set_title("Affinity vs expression (size ~ score)")
    cbar2 = plt.colorbar(sc, ax=ax2, fraction=0.046, pad=0.01)
    cbar2.set_label("oriented log2FC", fontsize=10)

    to_label = pick_labels(gs_sub, args.label_top, args.label_metric)
    for i, row in gs_sub.iterrows():
        if i in to_label:
            ax2.annotate(row["gene"], (row["log2_aff_ratio"], row["log2FC_oriented"]),
                         xytext=(5,5), textcoords="offset points", fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.12)
    plt.savefig(f"{args.outdir}/fig_tf_panel_{tf}.pdf")
    plt.close(fig)
