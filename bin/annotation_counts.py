#!/usr/bin/env python3
"""Explicit new-dataset mode: sum annotated gene-symbol counts without a historical selection."""
import argparse,csv,json,re
from pathlib import Path
import pandas as pd

def build(samplesheet,files):
    samples=[r['sample_name'] for r in csv.DictReader(Path(samplesheet).open(),delimiter='\t') if r['assay']=='RNA']
    tables={};metadata=None
    for path in files:
        name=Path(path).name.removesuffix('.featureCounts.txt')
        table=pd.read_csv(path,sep='\t',comment='#',index_col=0)
        if table.index.has_duplicates or name in tables:raise ValueError('Duplicate count identities')
        annotation=table[['gene_name','Chr','Strand','Length']]
        if metadata is not None and not annotation.equals(metadata):raise ValueError('RNA annotation identities differ')
        metadata=annotation;tables[name]=table.iloc[:,-1]
    if set(samples)!=set(tables):raise ValueError('RNA count/sample identities differ')
    raw=pd.DataFrame({name:tables[name] for name in samples})
    if raw.isna().any().any() or (raw<0).any().any() or (raw%1!=0).any().any():raise ValueError('Expected nonnegative integer gene counts')
    names=metadata['gene_name'].astype('string').str.upper()
    keep=names.notna() & ~names.isin(['','.','NA']) & ~names.str.match(r'^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$',na=False)
    if names[keep].str.contains(';',regex=False).any():raise ValueError('Ambiguous multiple gene symbols require an explicit annotation mapping')
    if not keep.any():raise ValueError('No annotated symbols')
    groups=names[keep].groupby(names[keep],sort=True).groups
    rows=[];representatives=[];ledger=[]
    for symbol,members in groups.items():
        members=sorted(members);representative=members[0]
        rows.append(raw.loc[members].sum());representatives.append(representative)
        for member in members:ledger.append(dict(gene_symbol=symbol,representative_gene_id=representative,member_gene_id=member,policy='sum_all_annotated_ids_per_uppercase_symbol'))
    counts=pd.DataFrame(rows,index=representatives,columns=samples);counts.index.name='Geneid';counts.to_csv('annotation_counts.csv')
    pd.DataFrame(ledger).to_csv('annotation_groups.tsv',sep='\t',index=False)
    mapping=pd.DataFrame({'gene_name':metadata.loc[keep,'gene_name']});mapping.index.name='gene_id';mapping.to_csv('annotation_idmap.tsv',sep='\t')
    lengths=pd.DataFrame({'merged':pd.to_numeric(metadata.loc[keep,'Length'])});lengths.index.name='gene'
    if (lengths.merged<=0).any():raise ValueError('Nonpositive merged exon length')
    lengths.to_csv('annotation_lengths.tsv',sep='\t')
    expected=raw.loc[keep].sum();observed=counts.sum()
    if not expected.equals(observed):raise ValueError('Annotated read-count totals were not conserved')
    record=dict(mode='annotation_symbol_sum',historical_equivalence='NOT_CLAIMED',n_raw_gene_ids=len(raw),n_kept_gene_ids=int(keep.sum()),n_symbols=len(counts),n_collapsed_ids=int(keep.sum())-len(counts),sample_order=samples,annotated_count_totals={k:int(v) for k,v in expected.items()},discarded_unmapped_name_count_totals={k:int(v) for k,v in raw.loc[~keep].sum().items()},count_conservation='EXACT',mapping='featureCounts gene_name from the supplied GTF',lengths='featureCounts merged exon lengths; the existing TPM port averages lengths per uppercase symbol',representative_id='lexicographically first ID, used only to index the summed symbol count')
    Path('annotation_counts_provenance.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--counts',nargs='+',required=True);a=p.parse_args();build(a.samplesheet,a.counts)
