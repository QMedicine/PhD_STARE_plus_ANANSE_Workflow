#!/usr/bin/env python3
"""Standalone port of the archived tfactivity mean-aggregation/TPM/filter steps."""
import argparse,json,platform
from pathlib import Path
import pandas as pd

def calculate(counts_file,lengths_file,map_file,min_count=50,min_tpm=0.1,min_count_tf=50,min_tpm_tf=0.1):
    mapping=pd.read_csv(map_file,sep='\t',index_col=0)
    lengths=pd.read_csv(lengths_file,sep='\t',index_col=0,usecols=['gene','merged'])
    lengths.columns=['length']
    raw=pd.read_csv(counts_file,index_col=0)
    mapping.index=mapping.index.map(lambda s:s.split('.')[0]);raw.index=raw.index.map(lambda s:s.split('.')[0])
    aliases=mapping['gene_name'].to_dict()
    existing=mapping['gene_name'].str.upper().tolist()
    index=raw.index.map(lambda x:aliases.get(x,x)).str.upper()
    missing=(~index.isin(existing)).sum()
    if len(raw) and missing/len(raw)>0.1:raise ValueError('More than 10% of genes absent from GTF mapping')
    raw.index=index;raw=raw.loc[raw.index.isin(existing)].apply(pd.to_numeric).groupby(level=0).mean()
    raw.to_csv('counts.clean.tsv',sep='\t')
    # Match the archived process boundary: calculate_tpm.py reopens the clean TSV.
    raw=pd.read_csv('counts.clean.tsv',sep='\t',index_col=0)
    lengths.index=lengths.index.map(aliases).str.upper();lengths=lengths.div(1000).groupby(level=0).mean()
    shared=lengths.index.intersection(raw.index);lengths=lengths.loc[shared];counts=raw.loc[shared]
    if (lengths['length']<=0).any():raise ValueError('Nonpositive gene lengths')
    rpk=counts.div(lengths['length'],axis=0);scales=rpk.sum()/1e6
    if (scales<=0).any():raise ValueError('Zero TPM library denominator')
    tpm=rpk.div(scales,axis=1);tpm.to_csv('counts.tpm.tsv',sep='\t')
    for prefix,mc,mt in [('counts',min_count,min_tpm),('TFs',min_count_tf,min_tpm_tf)]:
        selected=raw.loc[raw.sum(axis=1)>=mc].index.intersection(tpm.loc[tpm.mean(axis=1)>=mt].index)
        cf=raw.loc[selected].copy();tf=tpm.loc[selected].copy();cf.index.name=tf.index.name='gene_id'
        cf.to_csv(prefix+'.counts_filtered.tsv',sep='\t');tf.to_csv(prefix+'.tpm_filtered.tsv',sep='\t')
        Path(prefix+'.genes_filtered.txt').write_text('\n'.join(selected))
    Path('versions.json').write_text(json.dumps(dict(python=platform.python_version(),pandas=pd.__version__),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['counts','lengths','idmap']:p.add_argument('--'+key,required=True)
    for key,default in [('min-count',50),('min-count-tf',50),('min-tpm',0.1),('min-tpm-tf',0.1)]:p.add_argument('--'+key,type=float,default=default)
    a=p.parse_args();calculate(a.counts,a.lengths,a.idmap,a.min_count,a.min_tpm,a.min_count_tf,a.min_tpm_tf)
