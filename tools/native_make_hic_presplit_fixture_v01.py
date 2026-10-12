#!/usr/bin/env python3
"""Small deterministic split-read fixture; no historical-equivalence claim."""
import csv,gzip,json,random,subprocess
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/inputs/hic_presplit_fixture_v01')
root.mkdir(parents=True,exist_ok=True)
rng=random.Random(731);sequence=list(''.join(rng.choice('ACGT') for _ in range(30000)))
truth=[];reads=[[],[]]
for i in range(8):
    name=f'synthetic_hic_{i:02d}'
    for mate in [1,2]:
        pos=1000+i*1200+(300 if mate==2 else 0);other=15000+i*600+(200 if mate==2 else 0)
        split=(i%2==0 if mate==1 else i%3==0)
        if split:
            sequence[pos+40:pos+44]=list('GATC');sequence[other:other+4]=list('GATC')
        truth.append(dict(qname=name,mate=mate,position=pos,other=other,split=split))
sequence=''.join(sequence)
(root/'reference.fa').write_text('>test_contig\n'+sequence+'\n')
for row in truth:
    pos=row['position'];other=row['other']
    seq=sequence[pos:pos+44]+sequence[other:other+44] if row['split'] else sequence[pos:pos+88]
    if ('GATCGATC' in seq)!=row['split']:raise ValueError('Unplanned ligation signature')
    reads[row['mate']-1].append('@'+row['qname']+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
for mate in [1,2]:
    with gzip.open(root/f'fixture_R{mate}.fastq.gz','wt') as f:f.write(''.join(reads[mate-1]))
(root/'truth.json').write_text(json.dumps(truth,indent=2)+'\n')
(root/'index').mkdir(exist_ok=True)
subprocess.run(['bowtie2-build',str(root/'reference.fa'),str(root/'index/genome')],check=True)
fields=['sample_name','replicate','species','assay','paired','cell_type','condition','cell_state','fastq_1','fastq_2','bam','peak']
rows=[]
for state in ['State_1','State_2']:
    for assay in ['RNA','ATAC','HiC']:
        for replicate in range(1,2 if assay=='HiC' else 3):
            row={k:'' for k in fields};row.update(sample_name=f'{state}_{assay}_{replicate}',replicate=str(replicate),species='synthetic',assay=assay,paired='TRUE' if assay=='HiC' else 'FALSE',cell_type='Synthetic_test',condition=state,cell_state=state,bam=str(root/'reference.fa'))
            if assay=='HiC':row.update(bam='',fastq_1=str(root/'fixture_R1.fastq.gz'),fastq_2=str(root/'fixture_R2.fastq.gz'))
            rows.append(row)
with (root/'samplesheet.tsv').open('w') as f:
    w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
(root/'source2target.tsv').write_text('source\ttarget\nState_1\tState_2\n')
task=root.parent.parent
params=dict(stage='raw_prepare',samplesheet=str(root/'samplesheet.tsv'),source2target=str(root/'source2target.tsv'),raw_sample_names=['State_1_HiC_1'],hic_alignment_mode='diffhic_presplit_standard_compat',hic_restriction_sites=['GATC'],hic_ligation_signature='GATCGATC',bowtie_index=str(root/'index'),bowtie_index_name='genome',raw_tools_container=str(task/'containers/raw_tools_v07.sif'),hic_raw_container=str(task/'containers/hic_raw_standard_v01.sif'),ananse_container=str(task/'containers/ananse_v03.sif'),container_binds=[str(task)])
(root/'params.json').write_text(json.dumps(params,indent=2)+'\n')
(root/'fixture_identity.json').write_text(json.dumps(dict(status='SYNTHETIC',pairs=8,split_segments=sum(r['split'] for r in truth),seed=731,placeholder_bam_note='Unused RNA/ATAC entries reference the synthetic FASTA solely for manifest file-existence validation; no assay biology is simulated'),indent=2)+'\n')
