#!/usr/bin/env python3
# 10_compute_gene_gabc_ratio.py  (robust header detection + --agg + optional overrides)
import argparse, re, pandas as pd, numpy as np

p = argparse.ArgumentParser()
p.add_argument("--abc_old", nargs="+", required=True)
p.add_argument("--abc_young", nargs="+", required=True)
p.add_argument("--agg", choices=["sum","mean","median"], default="mean")
p.add_argument("--out", required=True)
# Optional explicit columns if autodetect ever fails
p.add_argument("--gene-col", default=None, help="Exact column name for gene (e.g., 'Gene Name' or 'Ensembl ID')")
p.add_argument("--score-col", default=None, help="Exact column name for ABC score (e.g., 'ABC-Score')")
args = p.parse_args()

def norm(s: str) -> str:
    # normalize header keys: lower, remove spaces/_/-
    return re.sub(r"[ \t_\-]", "", s.lower())

GENE_PREFS = ["genename", "gene", "targetgene", "target_gene", "Gene Name"]   # normalized forms
GENE_FALLBACKS = ["ensemblid", "ensembl_id", "ensgid", "geneid", "Ensembl ID"]
SCORE_PREFS = ["abcscore", "abc_score", "gabc", "gabcscore", "gabc_score", "score", "ABC-Score"]

def pick_columns(df: pd.DataFrame, gene_override=None, score_override=None):
    lcmap = {norm(c): c for c in df.columns}
    if gene_override:
        gene_col = gene_override
        if gene_col not in df.columns:
            raise RuntimeError(f"Gene column '{gene_col}' not found. Available: {list(df.columns)}")
    else:
        gene_col = next((lcmap[k] for k in GENE_PREFS if k in lcmap), None)
        if gene_col is None:
            gene_col = next((lcmap[k] for k in GENE_FALLBACKS if k in lcmap), None)
    if gene_col is None:
        raise RuntimeError("Could not find a gene column (tried Gene Name and Ensembl IDs).")

    if score_override:
        score_col = score_override
        if score_col not in df.columns:
            raise RuntimeError(f"Score column '{score_col}' not found. Available: {list(df.columns)}")
    else:
        score_col = next((lcmap[k] for k in SCORE_PREFS if k in lcmap), None)
    if score_col is None:
        raise RuntimeError("Could not find an ABC score column (e.g., 'ABC-Score').")

    return gene_col, score_col

def read_interactions(path, gene_override=None, score_override=None):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    gene_col, score_col = pick_columns(df, gene_override, score_override)
    sub = df[[gene_col, score_col]].copy()
    sub.columns = ["gene", "score"]
    # Harmonize IDs/symbols: strip version, uppercase symbols
    sub["gene"] = sub["gene"].astype(str).str.split(".").str[0].str.upper()
    return sub

def gene_aggregate(file_list, how, gene_override=None, score_override=None):
    # Sum per gene within each replicate, then aggregate across replicates
    per_rep = []
    for f in file_list:
        g = (read_interactions(f, gene_override, score_override)
             .groupby("gene", as_index=True)["score"].sum()
             .rename(f))
        per_rep.append(g)
    M = pd.concat(per_rep, axis=1).fillna(0.0)
    if how == "sum":    agg = M.sum(axis=1)
    elif how == "mean": agg = M.mean(axis=1)
    else:               agg = M.median(axis=1)
    return agg.to_frame("cond_score")

old = gene_aggregate(args.abc_old,  args.agg, args.gene_col, args.score_col)
you = gene_aggregate(args.abc_young, args.agg, args.gene_col, args.score_col)

common = old.join(you, how="inner", lsuffix="_old", rsuffix="_young")
# Small pseudocount avoids log2(0)
common["log2_gabc_ratio"] = np.log2((common["cond_score_old"] + 1e-12) /
                                    (common["cond_score_young"] + 1e-12))
common.to_csv(args.out, sep="\t")
