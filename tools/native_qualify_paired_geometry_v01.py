#!/usr/bin/env python3
"""Inspect actual paired BAM flags, fragment lengths and ATAC peak geometry."""
import argparse,gzip,json,subprocess
from pathlib import Path

def qualify(run,fixture,assay,name,output):
    run=Path(run);identity=json.loads((Path(fixture)/'fixture_identity.json').read_text())
    folder=run/'results/raw'/assay/name
    bam=folder/(name+'.no_blacklist.bam' if assay=='ATAC' else name+'_Aligned.sortedByCoord.out.bam')
    command=['samtools','view','-f','66','-F','2308',str(bam)]
    lines=subprocess.check_output(command,text=True).splitlines()
    if len(lines)<0.9*identity['fragments']:raise ValueError('Too few mapped, proper paired read1 alignments')
    lengths=[abs(int(line.split('\t')[8])) for line in lines]
    if any(n<=0 or not 250<=n<=450 for n in lengths):raise ValueError('Fragment lengths differ from declared synthetic geometry')
    record=dict(status='PAIRED_GEOMETRY_PASS_SYNTHETIC_ONLY',assay=assay,run=str(run),input_fragments=identity['fragments'],proper_mapped_read1=len(lines),observed_fragment_length_min=min(lengths),observed_fragment_length_max=max(lengths),generic_states=['State_1','State_2'])
    if assay=='ATAC':
        peaks=list((folder/'peaks').glob('*_summits.bed'))
        if len(peaks)!=1 or not peaks[0].read_text().strip():raise ValueError('Synthetic paired ATAC produced no peaks')
        macs_commands=list((run/'work').glob('*/*/.command.sh'))
        commands=[p.read_text() for p in macs_commands if 'macs3 callpeak' in p.read_text()]
        if len(commands)!=1 or '-f BAMPE' not in commands[0] or '--shift -100' in commands[0]:raise ValueError('MACS3 did not use paired-end geometry')
        record['summit_count']=len(peaks[0].read_text().splitlines());record['macs_format']='BAMPE'
    else:
        counts=[]
        for p in sorted(folder.glob('*.trimmed_R*.fastq.gz')):
            with gzip.open(p,'rt') as f:counts.append(sum(1 for _ in f)//4)
        if counts!=[identity['fragments']]*2:raise ValueError('Trimmed mates are missing or unbalanced')
        record['trimmed_reads_per_mate']=counts
    Path(output).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['run','fixture','assay','name','output']:p.add_argument('--'+key,required=True)
    qualify(**vars(p.parse_args()))
