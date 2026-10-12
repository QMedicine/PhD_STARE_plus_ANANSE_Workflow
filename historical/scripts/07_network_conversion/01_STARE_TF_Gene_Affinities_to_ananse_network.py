#!/usr/bin/env python3
import argparse, gzip, os, sys, json
import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
import re

EM_DASH = "—"  # U+2014 em dash —

def read_table(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", dtype=str)

def load_stare_gene_tf(paths, gene_col_candidates=(
    "gene","Gene","gene_id","Geneid","Ensembl","Ensembl ID","EnsemblID",
    "ensembl_id","symbol","gene_name","Gene Name"
)):
    dfs = []
    for p in paths:
        df = read_table(p)
        # find a gene column; fall back to first column
        gene_col = next((c for c in gene_col_candidates if c in df.columns), df.columns[0])
        df = df.set_index(gene_col, drop=False)
        # remove non-TF cols if present
        # also drop STARE’s helper columns that are not TFs
        df = df.drop(columns=["NumPeaks", "AvgPeakDistance", "AvgPeakSize"], errors="ignore")
        non_tf = {gene_col, "gene_name","GeneName","symbol","chr","start","end","NumPeaks","AvgPeakDistance","AvgPeakSize"}
        tf_cols = [c for c in df.columns if c not in non_tf]
        # numeric, NA→0
        df = df[tf_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        dfs.append(df)

    # intersect genes/TFs across replicates
    genes = set(dfs[0].index)
    tfs   = set(dfs[0].columns)
    for d in dfs[1:]:
        genes &= set(d.index)
        tfs   &= set(d.columns)
    genes = sorted(genes)
    tfs   = sorted(tfs)

    # average across replicates
    mats = [d.loc[genes, tfs].to_numpy(dtype=float) for d in dfs]
    avg  = np.mean(mats, axis=0)
    out  = pd.DataFrame(avg, index=genes, columns=tfs)  # rows=genes, cols=TFs
    return out

def load_tpm(path, gene_id_cols=(
    "gene","gene_id","Gene","Geneid","symbol","gene_name",
    "ENSEMBL","EnsemblID","Ensembl ID"
)):
    df = pd.read_csv(path, sep="\t", dtype=str)
    gene_col = next((c for c in gene_id_cols if c in df.columns), df.columns[0])
    df = df.set_index(gene_col, drop=False)
    # take mean TPM across any sample columns
    num = df.drop(columns=[gene_col], errors="ignore").apply(pd.to_numeric, errors="coerce")
    tpm = num.mean(axis=1).fillna(0.0)
    tpm.name = "tpm"
    return tpm  # index = genes (IDs or symbols)

def _upper_key_index(series: pd.Series):
    """Return a Series keyed by UPPERCASE(index) for matching, values unchanged."""
    return pd.Series(series.values, index=[str(i).upper() for i in series.index])

# Strip pandas' duplicate-column suffixes like ".1", ".2" when matching TF symbols
_SUFFIX_RE = re.compile(r"\.\d+$")
def _base_symbol(name: str) -> str:
    return _SUFFIX_RE.sub("", str(name))

def load_id_symbol_map(path):
    m = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    # normalize to UPPER for robust matching
    m["gene_id_uc"] = m["gene_id"].str.upper()
    m["gene_name_uc"] = m["gene_name"].str.upper()
    id2sym = dict(zip(m["gene_id_uc"], m["gene_name_uc"]))
    sym2id = dict(zip(m["gene_name_uc"], m["gene_id_uc"]))
    return id2sym, sym2id

def build_gene_alias_tpm(tpm: pd.Series, id2sym: dict, sym2id: dict) -> pd.Series:
    """
    Return a TPM series indexed by UPPERCASE aliases so either Ensembl IDs
    or symbols can match:
      - keep original (UPPER) index
      - add Ensembl->symbol alias (if available)
      - add symbol->Ensembl alias (if available)
    """
    alias = {}
    for orig in tpm.index.astype(str):
        v = float(tpm.loc[orig])
        uc = orig.upper()
        alias[uc] = v
        if uc in id2sym:
            alias[id2sym[uc]] = v
        if uc in sym2id:
            alias[sym2id[uc]] = v
    return pd.Series(alias)

def make_network(stare_gene_tf: pd.DataFrame,
                 tpm_gene: pd.Series,
                 tpm_tf: pd.Series,
                 out_dir: str,
                 out_basename: str = "network",
                 pair_sep: str = EM_DASH,
                 min_binding: float = 0.0):
    os.makedirs(out_dir, exist_ok=True)

    # 1) Per-TF min-max scale gABC-informed binding to [0,1]
    bind = stare_gene_tf.copy()
    for tf in bind.columns:
        v = bind[tf].to_numpy(dtype=float)
        if np.max(v) > np.min(v):
            bind[tf] = minmax_scale(v, (0.0, 1.0))
        else:
            bind[tf] = 0.0

    # 2) Expression modulation (light, preserves binding as primary signal)
    #    Match TF symbols case-insensitively; targets matched as-is (best effort).
    tpm_tf_uc   = _upper_key_index(tpm_tf)
    tf_bases_uc = [_base_symbol(tf).upper() for tf in bind.columns]
    tf_expr     = tpm_tf_uc.reindex(tf_bases_uc).fillna(0.0).to_numpy(dtype=float)
    tf_expr_s   = minmax_scale(tf_expr, (0.0, 1.0)) if tf_expr.max() > 0 else np.zeros_like(tf_expr)
    tf_expr_s   = pd.Series(tf_expr_s, index=bind.columns)

    # targets (genes) – match case-insensitively (and via Ensembl<->symbol aliases if provided)
    tpm_gene_uc = tpm_gene if all(isinstance(i, str) and i.isupper() for i in tpm_gene.index) else _upper_key_index(tpm_gene)
    gene_uc = [str(g).upper() for g in bind.index]
    tgt_expr_raw = tpm_gene_uc.reindex(gene_uc).fillna(0.0).to_numpy(dtype=float)
    tgt_expr_s   = minmax_scale(tgt_expr_raw, (0.0, 1.0)) if tgt_expr_raw.max() > 0 else np.zeros_like(tgt_expr_raw)
    tgt_expr_s   = pd.Series(tgt_expr_s, index=bind.index)

    # 3) Build 2-column ANANSE-style network: tf_target, prob
    #    prob = binding_scaled * (0.5+0.5*tf_expr_s) * (0.5+0.5*tgt_expr_s)
    rows_tf_target = []
    rows_prob      = []

    # iterate genes first for cache locality
    for gene in bind.index:
        gt = float(tgt_expr_s.get(gene, 0.0))
        boost_target = 0.5 + 0.5 * gt
        b_row = bind.loc[gene].to_dict()
        for tf, b in b_row.items():
            b = float(b)
            if b <= min_binding:
                continue
            te = float(tf_expr_s.get(tf, 0.0))
            prob = b * (0.5 + 0.5 * te) * boost_target
            rows_tf_target.append(f"{tf}{pair_sep}{gene}")
            rows_prob.append(prob)

    out = pd.DataFrame({"tf_target": rows_tf_target, "prob": rows_prob})
    out.sort_values("tf_target", inplace=True)

    # 4) Write outputs expected by your ANANSE version
    tsv_path = os.path.join(out_dir, f"{out_basename}.tsv")
    net_path = os.path.join(out_dir, f"{out_basename}.network")
    out.to_csv(tsv_path, sep="\t", index=False)
    out.to_csv(net_path, sep="\t", index=False)

    # metadata for provenance
    meta = {
        "note": "Derived from STARE Gene_TF matrices (gABC-informed), per-TF min-max, expression-modulated",
        "n_edges": int(out.shape[0]),
        "n_genes": int(stare_gene_tf.shape[0]),
        "n_tfs": int(stare_gene_tf.shape[1]),
        "columns": list(out.columns),
        "pair_separator": pair_sep,
        "min_binding_threshold": min_binding
    }
    with open(os.path.join(out_dir, "network_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"[OK] Wrote:\n  {tsv_path}\n  {net_path}\nEdges: {out.shape[0]}")
    return net_path

def main():
    ap = argparse.ArgumentParser(
        description="Convert STARE TF_Gene_Affinities (gABC) into ANANSE-style network with columns: tf_target, prob"
    )
    ap.add_argument("--stare", nargs="+", required=True,
                    help="One or more STARE Gene_TF_matrices files (txt/tsv[.gz]); replicates will be averaged.")
    ap.add_argument("--tpm", required=True,
                    help="TPM expression table (genes as rows; columns are samples—mean TPM is used).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--basename", default="network", help="Base name for outputs (default: network).")
    ap.add_argument("--sep", default=EM_DASH, help="Separator between TF and target (default: em dash '—').")
    ap.add_argument("--min-binding", type=float, default=0.0,
                    help="Discard edges with scaled binding ≤ this value (default: 0.0).")
    ap.add_argument("--idmap", default=None,
                help="Optional 2-column mapping (gene_id\\tgene_name) for Ensembl<->symbol aliasing during expression matching.")
    args = ap.parse_args()

    stare = load_stare_gene_tf(args.stare)
    tpm   = load_tpm(args.tpm)
    # Build a case-insensitive TPM vector for target genes with Ensembl<->symbol aliases
    if args.idmap:
        id2sym, sym2id = load_id_symbol_map(args.idmap)
        tpm_gene_for_targets = build_gene_alias_tpm(tpm, id2sym, sym2id)
    else:
        tpm_gene_for_targets = _upper_key_index(tpm)


    # Reuse same TPM vector for TFs and targets; TF matching is case-insensitive.
    make_network(stare, tpm_gene=tpm_gene_for_targets, tpm_tf=tpm,
             out_dir=args.outdir, out_basename=args.basename,
             pair_sep=args.sep, min_binding=args.min_binding)

if __name__ == "__main__":
    main()
