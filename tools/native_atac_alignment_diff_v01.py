#!/usr/bin/env python3
"""Compare stored alignments by coordinate without retaining full BAMs in memory."""
from collections import Counter,defaultdict
from itertools import groupby
from pathlib import Path
import hashlib,json,pysam
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
name='SRR13186609_Young_HSC_B1_Rep1'
paths=[Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered')/(name+'_filtered.bam'),root/'runs/raw_ATAC_pilot_v17/results/raw/ATAC'/name/(name+'.filtered.bam')]
changes=[defaultdict(list),defaultdict(list)];counts=[0,0]
def coordinate(r):return (r.reference_id,r.reference_start)
def payload(r):return (r.query_name,r.flag,r.mapping_quality,r.cigarstring,r.query_sequence,r.qual,tuple(sorted((k,str(v)) for k,v in r.tags if k!='PG')))
with pysam.AlignmentFile(paths[0],'rb') as a,pysam.AlignmentFile(paths[1],'rb') as b:
    if a.references!=b.references or a.lengths!=b.lengths:raise ValueError('Reference dictionaries differ')
    groups=[iter(groupby(a,coordinate)),iter(groupby(b,coordinate))];current=[next(g,None) for g in groups]
    while any(x is not None for x in current):
        pos=min(x[0] for x in current if x is not None);values=[]
        for i in [0,1]:
            if current[i] is not None and current[i][0]==pos:
                v=Counter(payload(r) for r in current[i][1]);counts[i]+=sum(v.values());current[i]=next(groups[i],None)
            else:v=Counter()
            values.append(v)
        for i in [0,1]:
            for p,n in (values[i]-values[1-i]).items():
                q,flag,mapq,cigar,sequence,quality,tags=p
                # Normalize to original sequenced orientation before hashing.
                if flag&16:sequence=sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1];quality=quality[::-1]
                changes[i][q].append(dict(contig=a.references[pos[0]],position0=pos[1],flag=flag,mapq=mapq,cigar=cigar,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),quality_sha256=hashlib.sha256(quality.encode()).hexdigest(),tags=tags,multiplicity=n))
names=set(changes[0])|set(changes[1]);both=set(changes[0])&set(changes[1])
def signature(row):return row['sequence_sha256'],row['quality_sha256']
same_input=sum(Counter(signature(r) for r in changes[0][q])==Counter(signature(r) for r in changes[1][q]) for q in both)
out=root/'evidence/atac_alignment_diff_v01';out.mkdir(parents=True,exist_ok=True)
receipt=dict(status='DIAGNOSTIC_NOT_EQUIVALENCE',paths=[str(p) for p in paths],alignment_counts=counts,differing_read_names=len(names),different_in_both=len(both),only_historical=len(set(changes[0])-set(changes[1])),only_rebuilt=len(set(changes[1])-set(changes[0])),both_with_identical_sequence_and_quality=same_input,changes=[dict(historical=changes[0].get(q,[]),rebuilt=changes[1].get(q,[]),qname=q) for q in sorted(names)])
(out/'alignment_diff.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='changes'},indent=2))
