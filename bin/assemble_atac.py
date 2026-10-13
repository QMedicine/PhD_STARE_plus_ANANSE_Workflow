#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path

def assemble(samplesheet,files,mode,output,header_order=None,data_order=None):
    with open(samplesheet) as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['assay']=='ATAC']
    # Explicit file/sample association; never use filesystem enumeration order.
    available={}
    for file in files:
        name=Path(file).name
        if not name.endswith('.window_counts.tsv'):raise ValueError('Unexpected ATAC count filename '+name)
        sample=name.removesuffix('.window_counts.tsv')
        if sample in available:raise ValueError('Duplicate ATAC count identity '+sample)
        available[sample]=file
    manifest_names=[r['sample_name'] for r in rows]
    if not manifest_names or len(set(manifest_names))!=len(manifest_names):raise ValueError('Duplicate or empty ATAC manifest identities')
    if set(manifest_names)&{'#chrom','start','end','peak_id'}:raise ValueError('ATAC sample identity collides with activity annotation columns')
    if set(available)!={r['sample_name'] for r in rows}:raise ValueError('ATAC counts and manifest differ')
    labels=[r['sample_name'] for r in rows]
    actual=list(labels)
    if mode=='historical_column_order':
        # The archived mismatch is described by explicit identities in the run configuration.
        for name,order in [('header_order',header_order),('data_order',data_order)]:
            if not order or len(order)!=len(labels) or set(order)!=set(labels):
                raise ValueError(name+' must contain every ATAC sample exactly once')
        labels=list(header_order);actual=list(data_order)
    elif mode!='corrected_column_order':raise ValueError('Invalid ATAC mode')
    data=[];identities=None
    for sample in actual:
        with open(available[sample]) as f:table=list(csv.reader(f,delimiter='\t'))
        if not table or any(len(r)!=5 for r in table):raise ValueError('ATAC counts require BED4 plus exactly one count column')
        keys=[r[:4] for r in table]
        if len({r[3] for r in keys})!=len(keys):raise ValueError('Duplicate peak identities')
        for row in table:
            if not row[-1].isdigit() or not row[1].isdigit() or not row[2].isdigit() or int(row[2])<=int(row[1]):raise ValueError('Invalid peak coordinates or nonnegative integer count')
        if identities is not None and identities!=keys:raise ValueError('Peak identities/order differ between samples')
        identities=keys;data.append([r[-1] for r in table])
    with open(output,'w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['peak_id']+labels)
        for i,row in enumerate(identities):w.writerow([row[3]]+[col[i] for col in data])
    with open('union_base.tsv','w') as f:
        csv.writer(f,delimiter='\t',lineterminator='\n').writerows(identities)
    Path('column_assignment.json').write_text(json.dumps(dict(mode=mode,header_samples=labels,data_samples=actual,columns_match_samples=labels==actual),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--counts',nargs='+',required=True);p.add_argument('--mode',required=True);p.add_argument('--output',required=True);p.add_argument('--header-order',default='');p.add_argument('--data-order',default='');a=p.parse_args();assemble(a.samplesheet,a.counts,a.mode,a.output,a.header_order.split(',') if a.header_order else None,a.data_order.split(',') if a.data_order else None)
