#!/usr/bin/env python3
# 11_summarize_top_tfs.py — ID-first, single-ID mapping; adds TPM & TF DE; robust gABC
import argparse, re
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--tf_rank", required=True)   # ...tf_ranking.annot.tsv (from step 9) index = tf column name; has 'dcg'
p.add_argument("--score",   required=True)   # ...ATAC-seq.score.tsv (genes x TFs)
p.add_argument("--coef",    required=True)   # ...ATAC-seq.filtered.tsv (index TF, col 'value')
p.add_argument("--ratio",   required=True)   # ...ATAC-seq.tsv (B/A ratio: genes x TFs)
p.add_argument("--deseq",   required=True)   # DESeq2 results (index gene; has log2FoldChange, padj, optionally gene_name)
p.add_argument("--gabc",    required=False)  # ...log2_gabc_ratio.tsv (index = gene SYMBOL (UPPER), col 'log2_gabc_ratio')
p.add_argument("--delta",   type=float, default=0.25)
p.add_argument("--topN",    type=int,   default=20)
p.add_argument("--topM",    type=int,   default=25)
p.add_argument("--tf_tpm",  required=True)   # TFs.tpm_filtered.tsv (first col = gene_id (symbol) + replicate TPMs)
p.add_argument("--idmap",   default=None)    # id_symbol_map.txt (gene_id  gene_name)
p.add_argument("--expanded-map", default=None)  # optional jaspar_expanded_map.csv
p.add_argument("--tf_map",  required=True)   # tf_mapping.tsv from step 9 (one row per TF column)
p.add_argument("--out",     required=True)
args = p.parse_args()

# -------------------------------
# Load core inputs
# -------------------------------
tf_rank = pd.read_csv(args.tf_rank, sep="\t", index_col=0)     # index = tf column name; has 'dcg'
scores  = pd.read_csv(args.score,   sep="\t", index_col=0)     # genes x TFs (score)
coefs   = pd.read_csv(args.coef,    sep="\t", index_col=0)     # index TF, col 'value'
ratio   = pd.read_csv(args.ratio,   sep="\t", index_col=0)     # genes x TFs (B/A ratio)
de      = pd.read_csv(args.deseq,   sep="\t", index_col=0) # set index as 1st column OR
#de      = pd.read_csv(args.deseq,   sep="\t")
#de = de.set_index("gene_name") # Set the index of 'de' to be the "gene_name" column

# Strip Ensembl version in DE index if present
de.index = de.index.astype(str).str.split(".").str[0]

# Optional ID map to build symbol<->Ensembl helpers
sym2ens, ens2sym = {}, {}
if args.idmap:
    idmap = pd.read_csv(args.idmap, sep="\t", index_col=0)
    # Set the index of 'idmap' to be the "Geneid" column
    # idmap = idmap.set_index("Geneid")
    idmap.index = idmap.index.astype(str).str.split(".").str[0]
    idmap["gene_key"] = idmap["gene_name"].astype(str).str.upper()
    # symbol -> ensembl
    sym2ens = idmap.reset_index().set_index("gene_key").iloc[:,0].to_dict()
    # ensembl -> SYMBOL_UPPER
    ens2sym = idmap["gene_key"].to_dict()

# Align gene sets between scores and DE
genes = scores.index.intersection(de.index)
scores = scores.loc[genes]
de     = de.loc[genes]

# Optional gABC
gabc = None
if args.gabc:
    gabc = pd.read_csv(args.gabc, sep="\t", index_col=0)
    gabc.index = gabc.index.astype(str).str.upper()

# TF mapping (source of truth for TF ID and names)
# Expected cols (from step 9): tf_col, tf_base, jaspar_id, tf_symbol_clean, dimer, monomer, from_id, width
tfmap = pd.read_csv(args.tf_map, sep="\t")

# -------------------------------
# Pick TFs: topN by DCG and present in all matrices
# -------------------------------
rank_sorted = tf_rank.sort_values("dcg", ascending=False)
tfs = []
for tf in rank_sorted.index.tolist():
    if tf in scores.columns and tf in ratio.columns and tf in coefs.index:
        tfs.append(tf)
    if len(tfs) >= args.topN:
        break

# -------------------------------
# Build (tf, gene) records
# -------------------------------
rows = []
for tf in tfs:
    coef_val  = float(coefs.at[tf, "value"]) if tf in coefs.index else np.nan
    top_genes = scores[tf].nlargest(args.topM).index
    for g in top_genes:
        log2fc = float(de.at[g, "log2FoldChange"]) if "log2FoldChange" in de.columns else np.nan
        # convert B/A ratio to log2
        aff = float(ratio.at[g, tf]) if (g in ratio.index and tf in ratio.columns) else np.nan
        log2_aff = np.log2(aff) if (np.isfinite(aff) and aff > 0) else np.nan

        # gABC by gene SYMBOL (map Ensembl -> symbol if needed)
        log2_gabc = np.nan
        if gabc is not None:
            sym_up = ens2sym.get(g, None)
            if sym_up and sym_up in gabc.index and "log2_gabc_ratio" in gabc.columns:
                log2_gabc = float(gabc.at[sym_up, "log2_gabc_ratio"])

        # driver type
        da = abs(log2_aff)  if np.isfinite(log2_aff)  else -np.inf
        dg = abs(log2_gabc) if np.isfinite(log2_gabc) else -np.inf
        if np.isfinite(da) and np.isfinite(dg):
            if   da > dg + args.delta: dt = "motif-dominant"
            elif dg > da + args.delta: dt = "contact-dominant"
            else:                      dt = "mixed"
        elif np.isfinite(da): dt = "motif-only"
        elif np.isfinite(dg): dt = "contact-only"
        else:                 dt = "NA"

        rows.append({
            "tf": tf,
            "gene": g,
            "coef_value": coef_val,
            "log2FC": log2fc,
            "score": float(scores.at[g, tf]),
            "log2_aff_ratio": log2_aff,
            "log2_gabc_ratio": log2_gabc,
            "consistency_tf_gene": (np.sign(coef_val) == np.sign(log2fc)) if (np.isfinite(coef_val) and np.isfinite(log2fc)) else np.nan,
            "affinity_higher_in": "Old" if (np.isfinite(log2_aff) and log2_aff > 0) else ("Young" if np.isfinite(log2_aff) and log2_aff < 0 else "NA"),
            "driver_type": dt
        })

out = pd.DataFrame(rows)

# -------------------------------
# Merge the per-column TF mapping (ID-first; single exact jaspar_id)
# -------------------------------
for col in ["jaspar_id", "tf_symbol"]:  # drop any stale columns if present
    if col in out.columns:
        out.drop(columns=[col], inplace=True)

out = out.merge(tfmap, left_on="tf", right_on="tf_col", how="left")

def strip_variant(s: str) -> str:
    return re.sub(r"\.\d+$", "", str(s))

# Display name (for plots/labels): prefer monomer > clean > base > tf
def best_symbol_for_display(row):
    for k in ("monomer", "tf_symbol_clean", "tf_base"):
        v = row.get(k)
        if pd.notnull(v) and str(v).strip():
            return str(v)
    return strip_variant(row["tf"])

# Lookup name for TPM/DE (must be a true gene symbol): prefer monomer > clean > base
def best_symbol_for_lookup(row):
    for k in ("monomer", "tf_symbol_clean", "tf_base"):
        v = row.get(k)
        if pd.notnull(v) and str(v).strip():
            return str(v)
    return strip_variant(row["tf"])

out["tf_symbol"]       = out.apply(best_symbol_for_display, axis=1)
out["tf_symbol_lookup"] = out.apply(best_symbol_for_lookup, axis=1)


# Choose display/lookup symbol (prefer monomer, else cleaned, else base, else tf)
def best_symbol(row):
    for k in ("monomer", "tf_symbol_clean", "tf_base", "tf"):
        v = row.get(k, None)
        if pd.notnull(v) and str(v).strip():
            return str(v)
    return str(row["tf"])

out["tf_symbol"] = out.apply(best_symbol, axis=1)
# 'jaspar_id' from tfmap is the exact PWM ID (may have .A/.B for dimers)

# -------------------------------
# Attach TPM per TF + TF DE stats (single merge; no overlapping pre-created cols)
# -------------------------------
tpm_df = pd.read_csv(args.tf_tpm, sep="\t")
first_col = tpm_df.columns[0]  # typically 'gene_id' (symbol)
tpm_df["__key__"] = tpm_df[first_col].astype(str).str.upper()
tpm_df = tpm_df.set_index("__key__", drop=False)

expr_cols = [c for c in tpm_df.columns if c not in ("__key__", first_col)]
def pick_cols(cols, pat):
    rx = re.compile(pat, re.IGNORECASE)
    return [c for c in cols if rx.search(c)]
old_cols   = pick_cols(expr_cols, r"\bold\b|\baged\b")
young_cols = pick_cols(expr_cols, r"\byoung\b|\byng\b")

def tf_de_lookup(sym_upper: str):
    # 1) gene_name column in DE (most common)
    if "gene_name" in de.columns:
        hits = de.index[de["gene_name"].astype(str).str.upper() == sym_upper]
        if len(hits) > 0:
            r = de.loc[hits[0]]
            return float(r.get("log2FoldChange", np.nan)), float(r.get("padj", np.nan))
    # 2) symbol in DE index
    if sym_upper in de.index:
        r = de.loc[sym_upper]
        return float(r.get("log2FoldChange", np.nan)), float(r.get("padj", np.nan))
    # 3) map symbol -> Ensembl then lookup
    ens = sym2ens.get(sym_upper, None)
    if isinstance(ens, str):
        ens_core = ens.split(".")[0]
        if ens_core in de.index:
            r = de.loc[ens_core]
            return float(r.get("log2FoldChange", np.nan)), float(r.get("padj", np.nan))
    return np.nan, np.nan

def fetch_tpm(sym_upper: str):
    key = sym_upper
    # your TPM table is keyed by SYMBOL by default; if not found, try Ensembl via idmap
    if key not in tpm_df.index and sym2ens:
        ens = sym2ens.get(sym_upper, "")
        key = ens.upper() if isinstance(ens, str) else ""
    if key not in tpm_df.index:
        return {}
    row = tpm_df.loc[key]
    d = {f"TFTPM|{c}": pd.to_numeric(row.get(c, np.nan), errors="coerce") for c in expr_cols}
    y_vals = pd.to_numeric(row[young_cols], errors="coerce").dropna().values if young_cols else np.array([])
    o_vals = pd.to_numeric(row[old_cols],   errors="coerce").dropna().values if old_cols   else np.array([])
    d["TFTPM_mean_YOUNG"] = float(np.mean(y_vals)) if y_vals.size else np.nan
    d["TFTPM_mean_OLD"]   = float(np.mean(o_vals)) if o_vals.size else np.nan
    return d

tfs_unique = out["tf"].unique().tolist()
tpm_records = []
for tfcol in tfs_unique:
    lookup_sym = str(out.loc[out["tf"] == tfcol, "tf_symbol_lookup"].iloc[0]).upper()
    tpm_dict = fetch_tpm(lookup_sym)         # use lookup_sym
    lfc, padj = tf_de_lookup(lookup_sym)     # use lookup_sym
    tpm_records.append({"tf": tfcol, **tpm_dict, "TF_log2FC": lfc, "TF_padj": padj})

tpm_per_tf = pd.DataFrame(tpm_records).set_index("tf")
out = out.merge(tpm_per_tf, left_on="tf", right_index=True, how="left")

# -------------------------------
# Optional: expanded-map context aggregated by symbol (em_*)
# -------------------------------
if args.expanded_map:
    em = pd.read_csv(args.expanded_map)
    em["tf_key"] = em["tf"].astype(str).str.upper()

    def agg_unique(series):
        vals = [str(x) for x in series if pd.notnull(x) and str(x) not in ("nan", "NA", "")]
        return ";".join(sorted(set(vals))) if vals else ""

    em_agg = (em.groupby("tf_key")
                .agg(em_id=("id", agg_unique),
                     em_dimer=("dimer", agg_unique),
                     em_monomer=("monomer", agg_unique),
                     em_from_id=("from_id", agg_unique),
                     em_width=("width", agg_unique))
                .reset_index())

    out["__tf_key__"] = out["tf_symbol"].astype(str).str.upper()
    out = out.merge(em_agg, left_on="__tf_key__", right_on="tf_key", how="left")
    out.drop(columns=["__tf_key__", "tf_key"], inplace=True, errors="ignore")

# -------------------------------
# Label for plotting: "TF [JASPAR_ID]"
# -------------------------------
def make_label(row):
    sym = str(row.get("tf_symbol") or row.get("tf") or "")
    jid = str(row.get("jaspar_id") or "")
    return f"{sym} [{jid}]" if jid else sym

out["tf_label"] = out.apply(make_label, axis=1)

# -------------------------------
# Write
# -------------------------------
out.to_csv(args.out, sep="\t", index=False)
