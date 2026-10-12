# 00_make_id_symbol_map.py
#!/usr/bin/env python3
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

pat_id = re.compile(r'gene_id "([^"]+)"')
pat_nm = re.compile(r'gene_name "([^"]+)"')
seen = {}
with open(args.gtf) as fh, open(args.out, "w") as out:
    out.write("gene_id\tgene_name\n")
    for line in fh:
        if line.startswith("#"): continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene": continue
        attrs = fields[8]
        m_id = pat_id.search(attrs)
        m_nm = pat_nm.search(attrs)
        if not m_id or not m_nm: continue
        gid = m_id.group(1)
        gnm = m_nm.group(1)
        if gid not in seen:
            out.write(f"{gid}\t{gnm}\n")
            seen[gid] = 1
