#!/usr/bin/env python3
"""Top-ranked TF TPM and directional DE, using the network's sample selection."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
from select_state_tpm import read_table, state_vector

def prepare(ranking, de, tpm, samplesheet, source, target, output, settings,
            tpm_limits=(0,10), fc_limits=(-6,6), font_file=None):
    if source == target:
        raise ValueError('Source and target must differ')
    for limits in (tpm_limits, fc_limits):
        if len(limits)!=2 or not all(math.isfinite(float(x)) for x in limits) or limits[0]>=limits[1]:
            raise ValueError('Heatmap limits require two increasing finite values')
    if fc_limits[0]>=0 or fc_limits[1]<=0:
        raise ValueError('logFC limits must span zero')
    _, ranks = read_table(ranking)
    ranks = ranks[:20]
    factors = [r['factor'] for r in ranks]
    if len(factors)!=20 or len(set(factors))!=20:
        raise ValueError('Exactly 20 distinct ranked TFs are required')
    _, _, _, src = state_vector(tpm, samplesheet, source)
    _, _, _, tgt = state_vector(tpm, samplesheet, target)
    def aliases(vector):
        mapped = {}
        for gene,value in vector.items():
            key=gene.upper()
            if key in mapped:
                raise ValueError('Ambiguous case-insensitive TPM identity: '+gene)
            mapped[key]=value
        return mapped
    src,tgt=aliases(src),aliases(tgt)
    fields, rows = read_table(de)
    if not {'gene_id','log2FoldChange'}.issubset(fields):
        raise ValueError('Directional DE requires gene_id and log2FoldChange')
    dec = {}
    for row in rows:
        key=row['gene_id'].upper()
        if key in dec:raise ValueError('Duplicate directional DE gene identifier')
        dec[key]=row
    table=[]
    for i,row in enumerate(ranks,1):
        tf=row['factor'];key=re.sub(r'\.\d+$','',tf).upper()
        fc=dec.get(key,{}).get('log2FoldChange','NA')
        record=dict(rank=i,factor=tf,source=source,target=target,
                    source_mean_TPM=src.get(key,'NA'),target_mean_TPM=tgt.get(key,'NA'),
                    log2FoldChange=fc,padj=dec.get(key,{}).get('padj','NA'),
                    influence_score=row['influence_score'],
                    target_score_scaled=row['target_score_scaled'],G_score_scaled=row['G_score_scaled'])
        for state_key,log_key in [('source_mean_TPM','source_log2_TPM_plus1'),('target_mean_TPM','target_log2_TPM_plus1')]:
            value=record[state_key]
            record[log_key]=math.log2(value+1) if value!='NA' else 'NA'
        table.append(record)
    with Path(output).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,list(table[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(table)
    Path(settings).write_text(json.dumps(dict(source=source,target=target,
        tpm_limits=list(tpm_limits),fc_limits=list(fc_limits),width_mm=90,height_mm=118,
        font_family='Arial',font_size_pt=6,font_file=font_file,
        n_tfs=20,clipping='Colour values clipped; table retains unmodified values',
        missing='NA in table; grey cells; no imputation as measured zero'),indent=2)+'\n')
    return table

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for key in ['ranking','de','tpm','samplesheet','source','target','output','settings']:
        parser.add_argument('--'+key,required=True)
    parser.add_argument('--tpm-limits',nargs=2,type=float,default=[0,10])
    parser.add_argument('--fc-limits',nargs=2,type=float,default=[-6,6])
    parser.add_argument('--font-file')
    prepare(**vars(parser.parse_args()))
