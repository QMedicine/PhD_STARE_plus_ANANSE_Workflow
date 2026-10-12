#!/usr/bin/env python3
"""Validate manifest identities and entry routes without guessing biology from paths."""
import argparse,csv,json,re
from pathlib import Path
REQUIRED=['sample_name','replicate','species','assay','paired','cell_type','condition','cell_state','fastq_1','fastq_2','bam','peak']
TOKEN=re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')

def read_tsv(path,required):
    with Path(path).open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):raise ValueError(f'{path}: required columns {required}')
        rows=list(reader)
    if not rows:raise ValueError(f'{path}: no rows')
    if any(None in r for r in rows):raise ValueError(f'{path}: extra fields')
    return rows

def validate(samplesheet,source2target,reverse=True,check_files=True,header_order=None,data_order=None,stage=None,process_fastq=False,raw_samples=None):
    samples=read_tsv(samplesheet,REQUIRED)
    identities=set();keys={};state_meta={};routes=[];pools={}
    for i,r in enumerate(samples,2):
        for key in REQUIRED[:8]:
            if not TOKEN.fullmatch(r[key] or ''):raise ValueError(f'Row {i}: invalid {key}')
        if r['sample_name'] in identities:raise ValueError('Duplicate sample_name '+r['sample_name'])
        identities.add(r['sample_name'])
        key=(r['cell_state'],r['assay'],r['replicate'])
        pool=r.get('pool_group','')
        if pool:
            if r['assay']!='HiC' or not TOKEN.fullmatch(pool):raise ValueError('pool_group is an optional safe Hi-C library-pooling identity')
            pool_meta=(r['cell_state'],r['replicate'],r['species'],r['cell_type'],r['condition'])
            if pool in pools and pools[pool]!=pool_meta:raise ValueError('Conflicting pool_group metadata')
            pools[pool]=pool_meta
        if key in keys and not (pool and keys[key]==pool and r['assay']=='HiC'):raise ValueError('Duplicate state/assay/replicate '+str(key))
        keys[key]=pool
        meta=(r['species'],r['cell_type'],r['condition'])
        if r['cell_state'] in state_meta and state_meta[r['cell_state']]!=meta:raise ValueError('Conflicting state metadata')
        state_meta[r['cell_state']]=meta
        if r['assay'] not in ['HiC','ATAC','RNA']:raise ValueError('assay must be HiC, ATAC or RNA')
        if r['paired'] not in ['TRUE','FALSE']:raise ValueError('paired must be TRUE or FALSE')
        if r['assay']=='HiC' and r['paired']!='TRUE':raise ValueError('Hi-C requires paired-end input')
        options=[name for name in ['fastq_1','bam','cooler'] if r.get(name)]
        if len(options)!=1:raise ValueError(f'{r["sample_name"]}: exactly one of fastq_1, bam, cooler required')
        route=options[0]
        if route=='cooler' and r['assay']!='HiC':raise ValueError('cooler is only valid for HiC')
        if r['assay']=='HiC' and stage in ['hic_diffhic_counts','atac_diffhic_hic_contacts','bam_diffhic_hic_contacts'] and route=='cooler':raise ValueError('diffHic routes require Hi-C BAM or explicitly enabled FASTQ entry')
        if route=='fastq_1' and (bool(r['fastq_2'])!=(r['paired']=='TRUE')):raise ValueError('fastq_2 and paired disagree')
        if route!='fastq_1' and r['fastq_2']:raise ValueError('fastq_2 without fastq_1')
        if r['peak'] and r['assay']!='ATAC':raise ValueError('peak is only valid for ATAC')
        for name in ['fastq_1','fastq_2','bam','peak','cooler']:
            value=r.get(name,'')
            if value:
                if not Path(value).is_absolute() or any(x in value for x in ['\n','\r','\t',"'",'"','`','$']):raise ValueError('Input paths must be safe absolute paths')
                if check_files and not Path(value).is_file():raise ValueError('Missing file '+value)
        if stage not in [None,'preflight'] and route=='fastq_1':
            if r['assay']=='HiC' and stage not in ['raw_prepare','hic_diffhic_counts','atac_diffhic_hic_contacts','bam_diffhic_hic_contacts']:raise ValueError('Raw Hi-C requires raw_prepare or an explicit diffHic route')
            if stage!='raw_prepare' and not process_fastq:raise ValueError('FASTQ entry requires --process_fastq or the raw_prepare pilot stage')
        routes.append(dict(sample_name=r['sample_name'],assay=r['assay'],entry=route))
    for name,order in [('historical_header_samples',header_order),('historical_data_samples',data_order)]:
        if order is not None:
            atac_names={r['sample_name'] for r in samples if r['assay']=='ATAC'}
            if len(order)!=len(atac_names) or set(order)!=atac_names:
                raise ValueError(name+' must contain each ATAC identity exactly once')
    transitions=[];seen=set();primary_transitions=[];primary_seen=set()
    for r in read_tsv(source2target,['source','target']):
        source,target=r['source'],r['target']
        if source==target or source not in state_meta or target not in state_meta:raise ValueError('Transition labels must match distinct samplesheet cell_state values exactly')
        if state_meta[source][0]!=state_meta[target][0]:raise ValueError('Source and target species differ')
        for state in [source,target]:
            assays={x['assay'] for x in samples if x['cell_state']==state}
            if assays!={'HiC','ATAC','RNA'}:raise ValueError('Each transition state requires HiC, ATAC and RNA')
            for assay in ['RNA','ATAC']:
                if sum(x['cell_state']==state and x['assay']==assay for x in samples)<2:raise ValueError('At least two biological replicates required for '+assay)
        primary_key=frozenset([source,target]) if reverse else (source,target)
        if primary_key not in primary_seen:
            primary_transitions.append(dict(source=source,target=target));primary_seen.add(primary_key)
        for pair in [(source,target)]+([(target,source)] if reverse else []):
            if pair not in seen:transitions.append(dict(source=pair[0],target=pair[1]));seen.add(pair)
    if raw_samples:
        eligible={r['sample_name'] for r in samples if r['fastq_1']}
        if stage!='raw_prepare' or len(raw_samples)!=len(set(raw_samples)) or not set(raw_samples).issubset(eligible):raise ValueError('raw_sample_names must identify unique FASTQ entries in raw_prepare')
    return dict(schema_version='0.1.0',samples=samples,transitions=transitions,primary_transitions=primary_transitions,routes=routes,files_checked=check_files)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--source2target',required=True);p.add_argument('--reverse-contrast',choices=['TRUE','FALSE'],default='TRUE');p.add_argument('--skip-file-check',action='store_true');p.add_argument('--output',required=True);p.add_argument('--header-order',default='');p.add_argument('--data-order',default='');p.add_argument('--stage');p.add_argument('--process-fastq',choices=['TRUE','FALSE'],default='FALSE')
    p.add_argument('--raw-samples',default='')
    a=p.parse_args();result=validate(a.samplesheet,a.source2target,a.reverse_contrast=='TRUE',not a.skip_file_check,a.header_order.split(',') if a.header_order else None,a.data_order.split(',') if a.data_order else None,a.stage,a.process_fastq=='TRUE',a.raw_samples.split(',') if a.raw_samples else None);Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
