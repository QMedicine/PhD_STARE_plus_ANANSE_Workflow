#!/usr/bin/env python3
import argparse,gzip
from pathlib import Path

def split_contacts(input_file,outdir,chromosomes=None):
    root=Path(outdir);root.mkdir()
    chromosomes=[c.removeprefix('chr') for c in chromosomes] if chromosomes else None
    handles={c:gzip.open(root/f'chr{c}.txt.gz','wt') for c in (chromosomes or [])};counts={c:0 for c in (chromosomes or [])};raw_names={}
    try:
        with open(input_file) as f:
            for line in f:
                parts=line.rstrip('\n').split('\t')
                if len(parts)!=7:raise ValueError('Expected seven GInteractions columns')
                a,b=parts[0].removeprefix('chr'),parts[3].removeprefix('chr')
                if a==b and (chromosomes is None or a in handles):
                    if a in raw_names and raw_names[a]!=parts[0]:raise ValueError('Mixed chromosome aliases in contact input')
                    raw_names[a]=parts[0]
                    if a not in handles:handles[a]=gzip.open(root/f'chr{a}.txt.gz','wt');counts[a]=0
                    handles[a].write('\t'.join([parts[1],parts[4],parts[6]])+'\n');counts[a]+=1
    finally:
        for h in handles.values():h.close()
    if not counts or any(n==0 for n in counts.values()):raise ValueError('Missing contacts on required chromosomes: '+str(counts))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--outdir',required=True);p.add_argument('--chromosomes',default='');a=p.parse_args();split_contacts(a.input,a.outdir,a.chromosomes.split(',') if a.chromosomes else None)
