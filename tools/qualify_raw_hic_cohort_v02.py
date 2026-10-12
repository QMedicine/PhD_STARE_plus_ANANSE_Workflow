#!/usr/bin/env python3
"""Per-library raw Hi-C integrity and primary/secondary BAM qualification."""
import argparse,csv,gzip,hashlib,json,re,subprocess
from pathlib import Path
import pysam
p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--task-index',type=int,required=True);p.add_argument('--raw-results',required=True);p.add_argument('--outdir',required=True);p.add_argument('--inputs-only',action='store_true');p.add_argument('--input-qualification-root');a=p.parse_args()
with open(a.samplesheet) as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['assay']=='HiC']
row=rows[a.task_index];out=Path(a.outdir)/row['sample_name'];out.mkdir(parents=True,exist_ok=True)
if a.input_qualification_root:
    receipt=Path(a.input_qualification_root)/row['sample_name']/'qualification.json'
    record=json.loads(receipt.read_text())
    if record['status']!='QUALIFIED_RAW_INPUTS' or record['sample']!=row:raise ValueError('Raw input qualification identities differ')
    for item in record['inputs']:
        path=Path(item['path'])
        if not path.is_file() or path.stat().st_size!=item['bytes']:raise ValueError('Qualified raw input changed size')
    n=record['paired_records']
    record['input_qualification_receipt']=str(receipt)
else:
    identity=[]
    for key in ['fastq_1','fastq_2']:
        path=Path(row[key]);digest=hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(8*1024*1024),b''):digest.update(chunk)
        subprocess.run(['gzip','-t',str(path)],check=True)
        identity.append(dict(mate=key,path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest(),gzip_crc='PASS'))
    n=0;lengths=[{},{}]
    with pysam.FastxFile(row['fastq_1']) as f1,pysam.FastxFile(row['fastq_2']) as f2:
        from itertools import zip_longest
        for r1,r2 in zip_longest(f1,f2):
            if r1 is None or r2 is None or re.sub('/[12]$','',r1.name)!=re.sub('/[12]$','',r2.name):raise ValueError('Hi-C mate record identities differ')
            for i,r in enumerate([r1,r2]):
                if r.quality is None or len(r.sequence)!=len(r.quality):raise ValueError('Missing/mismatched FASTQ qualities')
                length=str(len(r.sequence));lengths[i][length]=lengths[i].get(length,0)+1
            n+=1
    if n==0:raise ValueError('Empty library')
    record=dict(status='QUALIFIED_RAW_INPUTS',sample=row,paired_records=n,read_length_counts=lengths,inputs=identity,historical_t32_equivalence='UNKNOWN')
if not a.inputs_only:
    bam=Path(a.raw_results)/'raw/HiC'/row['sample_name']/(row['sample_name']+'.name_sorted.bam')
    subprocess.run(['samtools','quickcheck',str(bam)],check=True)
    primary=[0,0];mapped=[0,0];secondary=0;duplicates=0;bad_clips=0
    with pysam.AlignmentFile(bam,'rb') as f:
        if f.header.to_dict().get('HD',{}).get('SO')!='queryname':raise ValueError('BAM is not name-sorted')
        for r in f:
            if not r.is_paired or r.is_read1==r.is_read2 or r.is_supplementary:raise ValueError('Invalid Hi-C mate/segment flag')
            if r.is_secondary:
                secondary+=1
                if not r.is_unmapped and not any(op==5 for op,length in r.cigartuples):bad_clips+=1
            else:
                mate=0 if r.is_read1 else 1;primary[mate]+=1;mapped[mate]+=not r.is_unmapped;duplicates+=r.is_duplicate
    if primary!=[n,n] or bad_clips:raise ValueError('Primary read conservation or split hard clipping failed')
    digest=hashlib.sha256()
    with bam.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):digest.update(chunk)
    record.update(status='QUALIFIED_RAW_HIC_BAM_INTERNAL_CONTRACT',bam=str(bam),bam_sha256=digest.hexdigest(),primary_mates=primary,mapped_primary_mates=mapped,secondary_segments=secondary,duplicate_primary_records=duplicates,missing_secondary_hard_clips=bad_clips)
(out/'qualification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='inputs'},indent=2))
