#!/usr/bin/env python3
import argparse,gzip
from pathlib import Path

def split_contacts(input_file,outdir):
    root=Path(outdir);root.mkdir()
    chromosomes=[str(x) for x in range(1,20)]+['X']
    handles={c:gzip.open(root/f'chr{c}.txt.gz','wt') for c in chromosomes};counts={c:0 for c in chromosomes}
    try:
        with open(input_file) as f:
            for line in f:
                parts=line.rstrip('\n').split('\t')
                if len(parts)!=7:raise ValueError('Expected seven GInteractions columns')
                a,b=parts[0].removeprefix('chr'),parts[3].removeprefix('chr')
                if a==b and a in handles:
                    handles[a].write('\t'.join([parts[1],parts[4],parts[6]])+'\n');counts[a]+=1
    finally:
        for h in handles.values():h.close()
    if any(n==0 for n in counts.values()):raise ValueError('Missing contacts on required chromosomes: '+str(counts))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--outdir',required=True);a=p.parse_args();split_contacts(a.input,a.outdir)
