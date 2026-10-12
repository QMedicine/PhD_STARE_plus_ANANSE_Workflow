# 02_aggregate_synonyms.py
#!/usr/bin/env python3
import argparse, pandas as pd, re

parser = argparse.ArgumentParser()
parser.add_argument("--in", dest="inp", required=True)
parser.add_argument("--idmap", required=True)   # gene_id\tgene_name
parser.add_argument("--agg", choices=["mean","max","sum"], default="max")
parser.add_argument("--out", required=True)
args = parser.parse_args()

df = pd.read_csv(args.inp, sep="\t", index_col=0)

# Drop optional STARE summary columns if they exist
for col in ["NumPeaks","AvgPeakDistance","AvgPeakSize"]:
    if col in df.columns:
        df = df.drop(columns=[col])

# Map Ensembl IDs to symbols (uppercase), then aggregate duplicate symbols
idmap = pd.read_csv(args.idmap, sep="\t", index_col=0)
conv = idmap["gene_name"].to_dict()
df.index = df.index.map(lambda x: str(x).split(".")[0]).map(conv).str.upper()

df = getattr(df.groupby(df.index), args.agg)()

# Collapse TF synonyms: strip anything in parentheses, trim spaces, then aggregate across duplicate columns
newcols = df.columns.to_series().str.replace(r"\(.*\)", "", regex=True).str.strip()
df.columns = newcols
df = getattr(df.groupby(df.columns, axis=1), args.agg)()

df.to_csv(args.out, sep="\t")
