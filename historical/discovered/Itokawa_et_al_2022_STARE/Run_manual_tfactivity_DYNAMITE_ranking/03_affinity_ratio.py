# 03_affinity_ratio.py
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np

p = argparse.ArgumentParser()
p.add_argument("--a", required=True)
p.add_argument("--b", required=True)
p.add_argument("--out", required=True)
args = p.parse_args()

A = pd.read_csv(args.a, sep="\t", index_col=0)
B = pd.read_csv(args.b, sep="\t", index_col=0)

# Align by intersection of genes and TFs (as in tfactivity’s ratio step)
idx = A.index.intersection(B.index)
A = A.loc[idx]; B = B.loc[idx]
if not A.columns.equals(B.columns):
    # Ensure identical TF columns
    common = A.columns.intersection(B.columns)
    A, B = A[common], B[common]

# Avoid div-by-zero
A = A.replace(0, np.finfo(float).eps)

R = (B / A).replace([np.inf, -np.inf], np.nan).dropna()
R = R.round(5)
R.to_csv(args.out, sep="\t")
