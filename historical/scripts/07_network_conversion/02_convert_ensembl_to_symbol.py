#!/usr/bin/env python3
import argparse, gzip, io, os, re, sys
from collections import OrderedDict

EM_DASH = "—"  # U+2014 em dash
ID_VER_RE = re.compile(r'\.\d+$')  # strip possible version suffixes on IDs (rare)
TF_SUFFIX_RE = ID_VER_RE  # strip pandas-added TF suffixes like ".1", ".2"

def smart_open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def load_id_to_symbol_from_gtf(gtf_path):
    """
    Build dict: gene_id -> gene_name (UPPERCASE).
    If multiple entries exist, the first seen 'gene' feature wins.
    """
    id2sym = {}
    with smart_open(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            # GTF has 9 columns; attributes in col 9
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            attrs = parts[8]
            # Parse attributes: key "value"; pairs separated by ;
            # We want gene_id and gene_name
            # Handle both "gene_id" and "gene_name" robustly
            gid = None
            gname = None
            for field in attrs.split(";"):
                field = field.strip()
                if not field:
                    continue
                # expect: key "value"
                kv = field.split(" ", 1)
                if len(kv) != 2:
                    continue
                key = kv[0].strip()
                val = kv[1].strip().strip('"')
                if key == "gene_id":
                    gid = ID_VER_RE.sub("", val)  # strip version if present
                elif key == "gene_name":
                    gname = val
            if gid and gname and gid not in id2sym:
                id2sym[gid] = gname.upper()  # force uppercase
    return id2sym

def convert_file(in_path, out_path, id2sym, sep):
    """
    Read a 2-column file with header: tf_target <tab> prob
    Replace the right-hand gene token in tf_target with SYMBOL (UPPERCASE) when it looks like an Ensembl ID and is present in map.
    Also strip TF suffixes ('.1', '.2', ...) and collapse duplicate TF–target edges by keeping the max prob.
    Writes same header and order-of-first-occurrence to out_path.
    """
    n_in, n_conv, n_missing = 0, 0, 0

    # Use smart_open for both input and output so .gz is handled if present
    with smart_open(in_path, "rt") as fin:
        header = fin.readline()
        if not header:
            # nothing to do
            with smart_open(out_path, "wt") as fout:
                pass
            return (0, 0, 0)

        header_cols = header.rstrip("\n").split("\t")
        try:
            tf_target_idx = header_cols.index("tf_target")
            prob_idx = header_cols.index("prob")
        except ValueError:
            tf_target_idx, prob_idx = 0, 1

        # Collect transformed edges into an OrderedDict to keep first-seen order.
        # If a duplicate key appears, retain the maximum prob.
        rows = OrderedDict()

        for line in fin:
            if not line.strip():
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            tf_target = cols[tf_target_idx]
            prob_str = cols[prob_idx]

            # parse prob robustly
            try:
                prob_val = float(prob_str)
            except Exception:
                prob_val = 0.0

            if sep not in tf_target:
                # malformed; keep as-is but we still consider for max collapse on raw key
                key = tf_target
                prev = rows.get(key)
                rows[key] = max(prev, prob_val) if prev is not None else prob_val
                continue

            tf, gene_token = tf_target.split(sep, 1)

            # --- strip TF suffixes like ".1", ".2"
            tf = TF_SUFFIX_RE.sub("", tf)

            # --- map right-hand gene token to SYMBOL (UPPERCASE) when it's Ensembl; else uppercase
            gene_key = ID_VER_RE.sub("", gene_token)
            if gene_key.startswith("ENSMUSG") or gene_key.startswith("ENSG"):
                symbol = id2sym.get(gene_key)
                if symbol:
                    gene_sym = symbol  # already uppercased in loader
                    n_conv += 1
                else:
                    gene_sym = gene_token.upper()
                    n_missing += 1
            else:
                gene_sym = gene_token.upper()

            key = f"{tf}{sep}{gene_sym}"

            # collapse duplicates by keeping max prob
            prev = rows.get(key)
            rows[key] = max(prev, prob_val) if prev is not None else prob_val

    # --- write out after collapsing duplicates ---
    with smart_open(out_path, "wt") as fout:
        fout.write("tf_target\tprob\n")
        for key, p in rows.items():
            fout.write(f"{key}\t{p}\n")

    return (n_in, n_conv, n_missing)

def main():
    ap = argparse.ArgumentParser(
        description="Convert Ensembl IDs in ANANSE-style network files to UPPERCASE gene symbols using a GTF."
    )
    ap.add_argument("--gtf", required=True, help="GTF file with gene_id and gene_name attributes.")
    ap.add_argument("--sep", default=EM_DASH, help="Separator used in tf_target between TF and gene (default: em dash '—').")
    ap.add_argument("--inplace", action="store_true", help="Overwrite input files instead of writing *_SYMBOLS.*")
    ap.add_argument("--out-suffix", default="_SYMBOLS", help="Suffix for output files when not using --inplace.")
    ap.add_argument("inputs", nargs="+", help="One or more .network/.tsv files (optionally .gz).")
    args = ap.parse_args()

    id2sym = load_id_to_symbol_from_gtf(args.gtf)
    if not id2sym:
        print("ERROR: No (gene_id -> gene_name) mappings were parsed from the GTF.", file=sys.stderr)
        sys.exit(2)

    total = {"lines":0, "converted":0, "missing":0}
    for ip in args.inputs:
        base, ext = os.path.splitext(ip)
        # keep .gz correctly if present
        if ext == ".gz":
            base2, ext2 = os.path.splitext(base)  # .tsv.gz -> (file.tsv, .gz)
            outp = ip if args.inplace else f"{base2}{args.out_suffix}{ext2}.gz"
        else:
            outp = ip if args.inplace else f"{base}{args.out_suffix}{ext}"
        n_in, n_conv, n_missing = convert_file(ip, outp, id2sym, args.sep)
        total["lines"] += n_in
        total["converted"] += n_conv
        total["missing"] += n_missing
        print(f"[OK] {os.path.basename(ip)} -> {os.path.basename(outp)} | lines={n_in} converted={n_conv} missing_map={n_missing}")

    print(f"[SUMMARY] lines={total['lines']} converted={total['converted']} missing_map={total['missing']}")

if __name__ == "__main__":
    main()
