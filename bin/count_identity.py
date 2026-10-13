"""Validate identities before pandas can rename duplicate headings or align rows."""
import csv
from pathlib import Path

ANNOTATIONS=['Geneid','Chr','Start','End','Strand','Length','gene_name','gene_biotype']

def headings(path,delimiter):
    with Path(path).open(newline='') as stream:
        rows=csv.reader((line for line in stream if not line.startswith('#')),delimiter=delimiter)
        header=next(rows,[])
    if not header or any(not x for x in header) or len(set(header))!=len(header):
        raise ValueError('Empty or duplicate count-table headings: '+str(path))
    return header

def featurecounts(path,pd):
    header=headings(path,'\t')
    if header[:8]!=ANNOTATIONS or len(header)!=9:
        raise ValueError('Expected featureCounts annotation columns and exactly one BAM count column: '+str(path))
    table=pd.read_csv(path,sep='\t',comment='#',index_col=0)
    if table.index.has_duplicates or table.index.isna().any():raise ValueError('Duplicate or missing featureCounts gene identities')
    values=pd.to_numeric(table.iloc[:,-1],errors='raise')
    if values.isna().any() or (values<0).any() or (values%1!=0).any():raise ValueError('Expected nonnegative integer gene counts')
    return table,values
