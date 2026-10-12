#!/usr/bin/env python3
"""Check primary 5-prime coordinates, split flags, clipping, QNAMEs and mate fields."""
import json,hashlib
from pathlib import Path
import pysam
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
fixture=root/'inputs/hic_presplit_fixture_v01'
bam=root/'runs/hic_presplit_fixture_v18/results/raw/HiC/State_1_HiC_1/State_1_HiC_1.name_sorted.bam'
truth={(r['qname'],r['mate']):r for r in json.loads((fixture/'truth.json').read_text())}
primary={};secondary=[];records=[]
with pysam.AlignmentFile(bam,'rb') as f:
    for read in f:
        if read.is_unmapped or not read.is_paired or read.is_read1==read.is_read2:raise ValueError('Bad mapping or mate flags')
        key=(read.query_name,1 if read.is_read1 else 2)
        if key not in truth:raise ValueError('Unexpected read name/mate')
        row=truth[key]
        if read.is_secondary:
            if not row['split'] or read.reference_start!=row['other'] or read.is_supplementary:raise ValueError('Unexpected 3-prime alignment')
            secondary.append(key)
        else:
            if key in primary or read.reference_start!=row['position'] or read.is_supplementary:raise ValueError('5-prime primary coordinates differ')
            primary[key]=read
        if row['split'] and sum(n for op,n in read.cigartuples if op==5)!=44:raise ValueError('Split clipping differs')
        records.append(read)
if set(primary)!=set(truth) or len(secondary)!=sum(r['split'] for r in truth.values()):raise ValueError('Missing primary/split segments')
for key,read in primary.items():
    mate=primary[(key[0],3-key[1])]
    if read.mate_is_unmapped or read.next_reference_start!=mate.reference_start:raise ValueError('Mate information was not synchronized')
out=root/'evidence/hic_presplit_fixture_v01';out.mkdir(parents=True,exist_ok=True)
receipt=dict(status='PASS_SYNTHETIC_SPLIT_READ_GEOMETRY',pairs=8,primary_alignments=len(primary),secondary_alignments=len(secondary),total_alignments=len(records),all_primary_5prime_coordinates_exact=True,split_hard_clips_bp=44,mate_information_synchronized=True,bam_sha256=hashlib.sha256(bam.read_bytes()).hexdigest(),historical_t32_equivalence='UNKNOWN',counts_route='NOT_YET_QUALIFIED')
(out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
