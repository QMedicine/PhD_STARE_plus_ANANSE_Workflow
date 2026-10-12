# 04_dynamite_preprocess.py
#!/usr/bin/env python3
import argparse, pandas as pd
import sys # Import sys to allow exiting on error

p = argparse.ArgumentParser()
p.add_argument("--ratio", required=True)
p.add_argument("--deseq", required=True)
p.add_argument("--out", required=True)
args = p.parse_args()

# Load 'aff', using its 1st column (gene names) as the index
aff = pd.read_csv(args.ratio, sep="\t", index_col=0)

# Load 'de' WITHOUT setting an index column
#de  = pd.read_csv(args.deseq, sep="\t")

# Load 'de' by setting an index column as 1st column
de  = pd.read_csv(args.deseq, sep="\t", index_col=0)

# Check that the required "gene_name" column actually exists
#if "gene_name" not in de.columns:
#    print(f"Error: The --deseq file must contain a column named 'gene_name'.")
#    sys.exit(1)

# Set the index of 'de' to be the "gene_name" column
#de = de.set_index("gene_name")

# Clip gene version suffixes if using ensemble IDs for different transcripts or versions of the gene (ENSMUSG00000095041.2 etc)
# This now clips the gene names in aff.index and de.index
#aff.index = aff.index.to_series().str.split(".").str[0]
#de.index  = de.index.to_series().str.split(".").str[0]

# Intersection compares gene names to gene names
genes = aff.index.intersection(de.index)
aff = aff.loc[genes]
de  = de.loc[genes]

# Collapse duplicated IDs only averaging numeric columns
#aff = aff.groupby(aff.index).mean(numeric_only=True)
#de  = de.groupby(de.index).mean(numeric_only=True)

aff["Expression"] = (de["log2FoldChange"] > 0).astype(int)
aff.to_csv(args.out, sep="\t")