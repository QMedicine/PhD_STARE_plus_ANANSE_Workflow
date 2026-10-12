#!/usr/bin/env python3
"""Bounded header-format hypothesis test; never alter the original FASTQ/BAM."""
import json,pysam,subprocess
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
out=root/'evidence/atac_header_seed_probe_v01';out.mkdir(parents=True,exist_ok=True)
bam=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered/SRR13186609_Young_HSC_B1_Rep1_filtered.bam')
index=out/'historical.bai'
if not index.exists():subprocess.run(['samtools','index','-o',str(index),str(bam)],check=True)
diff=json.loads((root/'evidence/atac_alignment_diff_v01/alignment_diff.json').read_text())
chosen=[r for r in diff['changes'] if len(r['historical'])==1 and len(r['rebuilt'])==1][:20]
original={}
with pysam.AlignmentFile(bam,'rb',index_filename=str(index)) as f:
    for entry in chosen:
        a=entry['historical'][0]
        read=next(r for r in f.fetch(a['contig'],a['position0'],a['position0']+1) if r.query_name==entry['qname'])
        seq=read.get_forward_sequence();quality=read.qual[::-1] if read.is_reverse else read.qual
        original[read.query_name]=(seq,quality,read)
formats={
    'bare':lambda q,n:q,
    'spot_length':lambda q,n:q+' '+q.split('.')[-1]+' length='+str(n),
    'length':lambda q,n:q+' length='+str(n),
    'spot_mate':lambda q,n:q+' '+q.split('.')[-1]+'/1',
    'spot':lambda q,n:q+' '+q.split('.')[-1],
    'slash_mate':lambda q,n:q+'/1',
}
def result(read):return (read.reference_name,read.reference_start,read.flag,read.mapping_quality,read.cigarstring,tuple(sorted((k,str(v)) for k,v in read.tags if k!='PG')))
probe=[]
for name,format in formats.items():
    fq=out/(name+'.fastq')
    with fq.open('w') as f:
        for q,(seq,quality,read) in original.items():f.write('@'+format(q,len(seq))+'\n'+seq+'\n+\n'+quality+'\n')
    sam=out/(name+'.sam')
    with (out/(name+'.log')).open('w') as log:
        subprocess.run(['bowtie2','-p','2','--very-sensitive','-x','/vast/projects/Immune_Aging/Reference_genomes/mm39/Bowtie2/mm39','-U',str(fq),'-S',str(sam)],stderr=log,check=True)
    with pysam.AlignmentFile(sam,'r') as f:observed={r.query_name.removesuffix('/1'):r for r in f}
    matches=sum(q in observed and result(read)==result(observed[q]) for q,(_,_,read) in original.items())
    probe.append(dict(format=name,tested=len(original),exact_alignment_payload_matches=matches))
(out/'probe.json').write_text(json.dumps(dict(status='HEADER_FORMAT_HYPOTHESIS_TEST',comparisons=probe,note='A match is evidence for this header hypothesis on a bounded subset; it does not qualify the complete FASTQ route'),indent=2)+'\n');print(json.dumps(probe,indent=2))
