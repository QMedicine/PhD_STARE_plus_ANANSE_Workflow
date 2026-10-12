#!/usr/bin/env python3
"""Check matrix bin widths before assigning STARE's contact-resolution parameter."""
import argparse,json
from pathlib import Path
from hicmatrix import HiCMatrix

def validate(matrix,resolution):
    obj=HiCMatrix.hiCMatrix(matrix)
    widths=[int(end)-int(start) for chrom,start,end,*rest in obj.cut_intervals]
    if not widths or max(widths)!=resolution or any(w<=0 or w>resolution for w in widths):
        raise ValueError('Hi-C interval widths disagree with chosen_res')
    if resolution<=0 or any(int(start)%resolution!=0 for chrom,start,end,*rest in obj.cut_intervals):
        raise ValueError('Hi-C bins do not align to the requested resolution')
    Path('hic_resolution.json').write_text(json.dumps(dict(matrix=Path(matrix).name,chosen_res=resolution,n_bins=len(widths),max_bin_width=max(widths),status='PASS'),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--matrix',required=True);p.add_argument('--resolution',type=int,required=True);a=p.parse_args();validate(a.matrix,a.resolution)
