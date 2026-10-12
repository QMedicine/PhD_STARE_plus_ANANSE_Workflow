#!/usr/bin/env python3
# 09_build_tf_mapping_and_annotate.py  — consume-in-order, suffix-aware mapping
import argparse, re, gzip
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--inputs", nargs="+", required=True,
               help="Raw *_TF_Gene_Affinities_*.txt.gz from STARE")
p.add_argument("--motifs", required=True,
               help="motifs.converted.transfac used by STARE (filtered set, ORDER MATTERS)")
p.add_argument("--expanded_map_csv", default=None,
               help="Optional: jaspar_expanded_map.csv (to add dimer metadata)")
p.add_argument("--tf_ranking", required=True,
               help="...ATAC-seq.tf_ranking.tsv (index = TF column names)")
p.add_argument("--out_map", required=True, help="Output: tf_mapping.tsv (per TF column)")
p.add_argument("--out_tf_rank_annot", required=True, help="Output: annotated TF ranking")
args = p.parse_args()

NON_TF = {"NUMPEAKS","AVGPEAKDISTANCE","AVGPEAKSIZE","EXPRESSION"}

def strip_variant(s: str) -> str:
    return re.sub(r"\.\d+$", "", str(s))

def norm_sym(s: str) -> str:
    # match how headers are written (case-insensitive); keep dots (so X.1 stays X.1)
    return re.sub(r"[ _-]+", "", str(s).strip()).upper()

# 1) Gather actual TF column names as **pandas sees them** (so we capture .1/.2 suffixes)
tf_cols_present = []
for f in args.inputs:
    df0 = pd.read_csv(f, sep="\t", nrows=0)
    cols = df0.columns.tolist()[1:]  # drop geneID
    for c in cols:
        cu = norm_sym(c)
        if cu not in NON_TF and c not in tf_cols_present:
            tf_cols_present.append(c)
# Preserve order; pandas already deduped (TF, TF.1, TF.2, ...)

# 2) Parse TRANSFAC in file order → list of (id, tf_symbol_clean, base_upper)
motif_rows = []
cur_id, cur_name = None, None
with open(args.motifs) as fh:
    for line in fh:
        if line.startswith("ID "):
            cur_id = line.strip().split()[1]
            cur_name = None
        elif line.startswith("NA "):
            if cur_id is not None and cur_name is None:
                cur_name = line.strip()[3:].strip()
                motif_rows.append((cur_id, cur_name, norm_sym(cur_name)))
        elif line.startswith("//"):
            cur_id, cur_name = None, None

motifs_df = pd.DataFrame(motif_rows, columns=["id","tf_symbol_clean","tf_symbol_norm"])

# 3) Optional dimer metadata by ID (expanded map)
if args.expanded_map_csv:
    em = pd.read_csv(args.expanded_map_csv)
    em = em.rename(columns=str.lower)  # id, tf, dimer, monomer, from_id, width
    motifs_df = motifs_df.merge(
        em[["id","tf","dimer","monomer","from_id","width"]],
        on="id", how="left"
    )
else:
    for c in ("tf","dimer","monomer","from_id","width"):
        motifs_df[c] = pd.NA

# 4) Index motif rows *per base symbol* to consume in order
by_base = {}
for _, r in motifs_df.iterrows():
    by_base.setdefault(r["tf_symbol_norm"], []).append(r)
# Track how many we've consumed for each base symbol
used_idx = {k: 0 for k in by_base}

# 5) Build per-column mapping by consuming the next row for the base symbol
mapping_rows = []
for col in tf_cols_present:
    base = strip_variant(col)
    base_norm = norm_sym(base)
    rec = {
        "tf_col": col,              # exact column name (may carry .1/.2)
        "tf_base": base,            # base symbol without suffix
        "jaspar_id": pd.NA,
        "tf_symbol_clean": pd.NA,
        "tf_from_expanded": pd.NA,
        "dimer": pd.NA, "monomer": pd.NA, "from_id": pd.NA, "width": pd.NA
    }
    if base_norm in by_base and used_idx[base_norm] < len(by_base[base_norm]):
        r = by_base[base_norm][ used_idx[base_norm] ]
        used_idx[base_norm] += 1
        rec.update({
            "jaspar_id": r["id"],
            "tf_symbol_clean": r["tf_symbol_clean"],
            "tf_from_expanded": r.get("tf", pd.NA),
            "dimer": r.get("dimer", pd.NA),
            "monomer": r.get("monomer", pd.NA),
            "from_id": r.get("from_id", pd.NA),
            "width": r.get("width", pd.NA),
        })
    mapping_rows.append(rec)

map_df_out = pd.DataFrame(mapping_rows)
map_df_out.to_csv(args.out_map, sep="\t", index=False)

# 6) Annotate ranking by exact column name
rank = pd.read_csv(args.tf_ranking, sep="\t", index_col=0)
rank.index.name = "tf_col"
ann = rank.merge(map_df_out.set_index("tf_col"), left_index=True, right_index=True, how="left")

# label fallback if needed
ann["tf_symbol_clean"] = ann["tf_symbol_clean"].where(ann["tf_symbol_clean"].notna(), ann.index)
ann.to_csv(args.out_tf_rank_annot, sep="\t")
