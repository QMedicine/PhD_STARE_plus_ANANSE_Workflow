# 01_affinity_mean.py
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, sys

parser = argparse.ArgumentParser()
parser.add_argument("--method", choices=["mean","sum","ratio","rank"], default="mean")
parser.add_argument("--inputs", nargs="+", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

dfs = [pd.read_csv(p, sep="\t", index_col=0) for p in args.inputs]

# Keep intersection for mean, union+fill(0) for sum/rank (tfactivity behavior)
if args.method in ["sum","rank"]:
    idx = dfs[0].index
    cols = dfs[0].columns
    for df in dfs[1:]:
        idx = idx.union(df.index)
        cols = cols.union(df.columns)
    dfs = [df.reindex(idx).fillna(0).reindex(columns=cols).fillna(0) for df in dfs]
else:
    idx = dfs[0].index
    for df in dfs[1:]:
        idx = idx.intersection(df.index)
    dfs = [df.loc[idx] for df in dfs]

# Sanity
if not all(df.shape == dfs[0].shape for df in dfs):
    sys.exit(f"[ERR] Shape mismatch: {[df.shape for df in dfs]}")
if not all(dfs[0].index.equals(df.index) for df in dfs[1:]):
    sys.exit("[ERR] Row names mismatch.")
if not all(dfs[0].columns.equals(df.columns) for df in dfs[1:]):
    sys.exit("[ERR] Column names mismatch.")

if args.method == "mean":
    res = sum(dfs) / len(dfs)
elif args.method == "sum":
    res = sum(dfs)
elif args.method == "rank":
    res = 1 - (sum(dfs).rank(ascending=False) / len(dfs[0].index))
else:
    sys.exit("[ERR] 'ratio' not valid here; use 03_affinity_ratio.py")

res.to_csv(args.out, sep="\t", index=True)
