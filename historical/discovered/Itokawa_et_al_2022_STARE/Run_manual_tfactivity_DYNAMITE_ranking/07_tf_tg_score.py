# 07_tf_tg_score.py
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np

p = argparse.ArgumentParser()
p.add_argument("--deseq", required=True)
p.add_argument("--aff", required=True)    # ratio table from step 3
p.add_argument("--coefs", required=True)  # filtered DYNAMITE coeffs (TF,value)
p.add_argument("--out", required=True)
a = p.parse_args()

# Load 'de' WITHOUT setting an index column
#de  = pd.read_csv(a.deseq, sep="\t")
# Set the index of 'de' to be the "gene_name" column
#de = de.set_index("gene_name")

# Load 'de' WITH setting an index column to the 1st column
de  = pd.read_csv(a.deseq, sep="\t", index_col=0)

aff = pd.read_csv(a.aff,  sep="\t", index_col=0)
cof = pd.read_csv(a.coefs, sep="\t", index_col=0)

# Harmonize gene IDs and TF names
de.index  = de.index.to_series().str.split(".").str[0]
aff.index = aff.index.to_series().str.split(".").str[0]

genes = de.index.intersection(aff.index)
assert len(genes) > 0, "No overlapping genes between DE and affinities."

de  = de.loc[genes].groupby(level=0).mean(numeric_only=True)
aff = aff.loc[genes].groupby(level=0).mean(numeric_only=True)

tfs = aff.columns.intersection(cof.index)
assert len(tfs) > 0, "No overlapping TFs between affinities and coefficients."

aff = aff[tfs]
cof = cof.loc[tfs]

# Score = |log2FC| × affinity × |coef|
scores = aff.mul(de["log2FoldChange"].abs(), axis=0).mul(cof["value"].abs(), axis=1)

assert not scores.empty, "TF–TG score matrix is empty."
scores.to_csv(a.out, sep="\t")
