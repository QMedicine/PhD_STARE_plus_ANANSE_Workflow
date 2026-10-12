# 06_filter_coefficients.py
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np

p = argparse.ArgumentParser()
p.add_argument("--in", dest="inp", required=True)
p.add_argument("--thr", type=float, default=0.1)
p.add_argument("--out", required=True)
a = p.parse_args()

df = pd.read_csv(a.inp, sep="\t")
# Expect columns: TF, value  (as written by DYNAMITE.R)
df = df[(df["value"] >= a.thr) | (df["value"] <= -a.thr)]
df.to_csv(a.out, sep="\t", index=False)
