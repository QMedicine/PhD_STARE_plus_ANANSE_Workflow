#!/usr/bin/env python3
"""Make an audited compatibility adapter from the unchanged pinned diffHic source."""
import argparse,hashlib,json
from pathlib import Path

PATCHES = [
    ("mode='rU'", "mode='r'", 'Python text-mode compatibility'),
    ('read5.is_read2=read5.is_read2=not is_first', 'read3.is_read2=read5.is_read2=not is_first', 'Assign mate-2 flag to both split segments'),
    ('read3.is_supplementary=True', 'read3.is_secondary=True', 'Mark 3-prime segment secondary for the documented diffHic BAM contract'),
    ('+cut_err)', "+cut_err.decode('utf-8', errors='replace'))", 'Decode cutadapt error bytes'),
    ('+map_err)', "+map_err.decode('utf-8', errors='replace'))", 'Decode Bowtie2 error bytes'),
    ('if map_proc.wait():\n        map_out, map_err=map_proc.communicate()', 'map_out, map_err=map_proc.communicate()\n    if map_proc.returncode:', 'Drain stderr before checking unsplit-mapper exit status'),
]

def prepare(source, output, audit):
    source=Path(source);raw=source.read_bytes()
    identity=json.loads(source.with_name('source_identity.json').read_text())
    expected=identity['files']['inst/python/presplit_map.py']['sha256']
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Pinned upstream mapper checksum differs')
    text=raw.decode();changes=[]
    for before,after,reason in PATCHES:
        count=text.count(before)
        expected_count=2 if before=='+map_err)' else 1
        if count!=expected_count:raise ValueError('Unexpected upstream patch site: '+reason)
        text=text.replace(before,after);changes.append(dict(before=before,after=after,reason=reason,occurrences=count))
    Path(output).write_text(text)
    Path(audit).write_text(json.dumps(dict(method='diffhic_presplit_standard_compat',upstream_commit=identity['commit'],upstream_sha256=expected,adapter_sha256=hashlib.sha256(text.encode()).hexdigest(),changes=changes,historical_t32_equivalence='UNKNOWN'),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True);p.add_argument('--audit',required=True);a=p.parse_args();prepare(a.source,a.output,a.audit)
