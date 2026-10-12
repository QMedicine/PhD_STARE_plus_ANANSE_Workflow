#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import pandas as pd

def merge(selection,files,output):
    selected=pd.read_csv(selection,index_col=0)
    if selected.index.has_duplicates:raise ValueError('Historical selected Ensembl IDs are duplicated')
    tables={};metadata=None
    for file in files:
        sample=Path(file).name.removesuffix('.featureCounts.txt')
        table=pd.read_csv(file,sep='\t',comment='#',index_col=0)
        if table.index.has_duplicates:raise ValueError('featureCounts IDs are duplicated')
        annotations=table[['gene_name','Chr','Strand','Length']]
        if metadata is not None and not annotations.equals(metadata):raise ValueError('featureCounts annotation identities differ')
        metadata=annotations
        tables[sample]=table.iloc[:,-1]
    if set(tables)!=set(selected.columns):raise ValueError('Recount sample identities differ from historical selection')
    raw=pd.DataFrame({s:tables[s] for s in selected.columns})
    if not set(selected.index).issubset(raw.index):raise ValueError('Selected historical genes absent from featureCounts')
    names=metadata['gene_name'].where(metadata['gene_name'].notna(),pd.Series(metadata.index,index=metadata.index))
    groups=names.groupby(names,sort=False).groups
    counts=raw.loc[selected.index].copy();ledger=[]
    for gid in selected.index:
        members=list(groups[names.loc[gid]])
        if len(members)==1:continue
        chromosomes={x for value in metadata.loc[members,'Chr'].astype(str) for x in value.split(';') if x}
        strands={x for value in metadata.loc[members,'Strand'].astype(str) for x in value.split(';') if x}
        same_locus=len(chromosomes)==1 and len(strands)==1
        included=members if same_locus else [gid]
        counts.loc[gid]=raw.loc[included].sum()
        if len(members)>1:
            for member in included:ledger.append(dict(retained_gene_id=gid,member_gene_id=member,group_symbol=names.loc[gid],action='sum_same_chromosome_and_strand' if same_locus else 'retain_checkpoint_representative'))
    pd.DataFrame(ledger,columns=['retained_gene_id','member_gene_id','group_symbol','action']).to_csv('cleanup_groups.tsv',sep='\t',index=False)
    counts.index.name='Geneid';counts.to_csv(output)
    mismatch=(counts.to_numpy()!=selected.to_numpy()).sum()
    Path('recount_comparison.json').write_text(json.dumps(dict(n_genes=len(counts),n_samples=len(tables),n_different_values=int(mismatch),historical_values_identical=bool(mismatch==0),cleanup='same chromosome and strand summed; distinct loci retain explicit selected ID'),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--selection',required=True);p.add_argument('--counts',nargs='+',required=True);p.add_argument('--output',required=True);a=p.parse_args();merge(a.selection,a.counts,a.output)
