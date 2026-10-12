#!/usr/bin/env python3
"""Qualify all RNA FASTQ producers against saved sequence/count and qualified RNA checkpoints."""
import argparse,gzip,hashlib,json
from pathlib import Path
import pandas as pd

def digest(path,decoded=False):
    h=hashlib.sha256();n=0
    opener=gzip.open if decoded else open
    with opener(path,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block);n+=len(block)
    return h.hexdigest(),n

def qualify(run,reference_run,reference_counts,trimmed_dir,output):
    run=Path(run);reference_run=Path(reference_run)
    manifest=json.loads((run/'results/provenance/validated_inputs.json').read_text())
    samples=[r for r in manifest['samples'] if r['assay']=='RNA']
    if not samples or any(not r['fastq_1'] or r['paired']!='FALSE' for r in samples):raise ValueError('This historical regression requires all RNA samples through single-end FASTQs')
    traces=list((run/'results/provenance').glob('trace_*.tsv'))
    with traces[0].open() as f:tasks=list(__import__('csv').DictReader(f,delimiter='\t'))
    if len(traces)!=1 or any(t['status']!='COMPLETED' or t['exit']!='0' for t in tasks):raise ValueError('RNA producers are incomplete')
    expected=pd.read_csv(reference_counts,sep='\t',comment='#',index_col=0)
    expected.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in expected.columns]
    records=[]
    for row in samples:
        name=row['sample_name'];raw=run/'results/raw/RNA'/name
        ea,en=digest(Path(trimmed_dir)/(name+'_trimmed.fastq.gz'),True);oa,on=digest(raw/(name+'.trimmed_R1.fastq.gz'),True)
        observed=pd.read_csv(run/'results/rna/featurecounts'/(name+'.featureCounts.txt'),sep='\t',comment='#',index_col=0)
        if set(expected.index)!=set(observed.index):raise ValueError('Raw count gene identities differ')
        x=expected[name];y=observed.iloc[:,-1].loc[x.index]
        records.append(dict(sample=name,trimmed_sequence_bytes_identical=ea==oa,expected_decoded_sha256=ea,observed_decoded_sha256=oa,expected_decoded_bytes=en,observed_decoded_bytes=on,n_featurecounts_genes=len(x),different_gene_count_values=int((x!=y).sum()),counts_exact=bool((x==y).all())))
    downstream=[]
    files=['recount/selected_counts.csv','tpm/counts.tpm.tsv','tpm/counts.counts_filtered.tsv']
    for t in manifest['primary_transitions']:files.append('de/'+t['source']+'_to_'+t['target']+'/forward.deseq2.tsv')
    for relative in files:
        a,_=digest(reference_run/'results/rna'/relative);b,_=digest(run/'results/rna'/relative)
        downstream.append(dict(path=relative,expected_sha256=a,observed_sha256=b,byte_identical=a==b))
    reverse=[]
    for t in manifest['primary_transitions']:
        folder=run/'results/rna/de'/(t['source']+'_to_'+t['target'])
        a=pd.read_csv(folder/'forward.deseq2.tsv',sep='\t',index_col=0);b=pd.read_csv(folder/'reverse.deseq2.tsv',sep='\t',index_col=0)
        exact=list(a.index)==list(b.index) and bool((a.log2FoldChange==-b.log2FoldChange).all())
        for c in a.columns:
            if c!='log2FoldChange':exact=exact and bool(((a[c]==b[c])|(a[c].isna()&b[c].isna())).all())
        reverse.append(dict(source=t['source'],target=t['target'],sign_reversal_exact=exact))
    receipt=dict(run=str(run),producer_tasks=len(tasks),raw_samples=records,downstream=downstream,reverse=reverse,status='EXACT_SEQUENCE_COUNTS_TPM_PRIMARY_DE' if all(r['counts_exact'] and r['trimmed_sequence_bytes_identical'] for r in records) and all(r['byte_identical'] for r in downstream) and all(r['sign_reversal_exact'] for r in reverse) else 'DIFFERENT')
    Path(output).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
    if receipt['status']=='DIFFERENT':raise SystemExit('RNA raw integration regression differs')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['run','reference-run','reference-counts','trimmed-dir','output']:p.add_argument('--'+key,required=True)
    qualify(**vars(p.parse_args()))
