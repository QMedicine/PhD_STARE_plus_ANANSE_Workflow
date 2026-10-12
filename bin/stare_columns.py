#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--activity',required=True);p.add_argument('--state',required=True);a=p.parse_args()
with open(a.samplesheet) as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['assay']=='ATAC' and r['cell_state']==a.state]
with open(a.activity) as f:header=f.readline().rstrip('\n').split('\t')
if not rows:raise ValueError('No ATAC samples in state')
print(','.join(str(header.index(r['sample_name'])+1) for r in rows))
