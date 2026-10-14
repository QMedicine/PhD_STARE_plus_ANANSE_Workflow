#!/usr/bin/env python3
"""Qualify this downstream-only state-TPM run before copying final products."""
import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale
from qualify_selected_reports_v01 import qualify

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''):digest.update(block)
    return digest.hexdigest()

def run(results,params,input_receipt,destination=None):
    root=Path(results);settings=json.loads(Path(params).read_text())
    bound=json.loads(Path(input_receipt).read_text())
    if settings['atac_modes']!=['corrected_column_order'] or settings['tpm_expression_mode']!='state_specific':
        raise ValueError('Only corrected ATAC / state-specific TPM is authorized for this run')
    for row in bound['inputs']:
        if sha(row['path'])!=row['sha256']:raise ValueError('Bound input changed: '+row['path'])
    mode='corrected_column_order';base=root/mode
    networks=[];expression=[]
    with Path(settings['samplesheet']).open() as handle:samples=list(csv.DictReader(handle,delimiter='\t'))
    for state in bound['states']:
        network=base/'networks'/state/(state+'_SYMBOLS.tsv')
        selected=json.loads((base/'networks'/state/'tpm_selection.json').read_text())
        expected=[r['sample_name'] for r in samples if r['assay']=='RNA' and r['cell_state']==state]
        if selected['tpm_expression_mode']!='state_specific' or selected['cell_state']!=state or selected['selected_rna_samples']!=expected:
            raise ValueError('Incorrect state-specific RNA selection')
        n=0
        for chunk in pd.read_csv(network,sep='\t',chunksize=500000):
            if chunk.columns.tolist()!=['tf_target','prob'] or not np.isfinite(chunk['prob']).all() or not chunk['prob'].between(0,1).all():raise ValueError('Invalid network score')
            n+=len(chunk)
        if n==0:raise ValueError('Empty network')
        networks.append(dict(state=state,path=str(network),edges=n,sha256=sha(network)))
        expression.append(selected)
    with Path(settings['de_manifest']).open() as handle:contrasts=list(csv.DictReader(handle,delimiter='\t'))
    reports=[];ranks=[]
    for contrast in contrasts:
        source,target=contrast['source'],contrast['target']
        folder=base/(source+'_to_'+target)/'tables'
        ranking=folder/'influence.txt';diffnet=folder/'influence_diffnetwork.tsv'
        table=pd.read_csv(ranking,sep='\t')
        if len(table)<20 or table.factor.duplicated().any():raise ValueError('Invalid ranking identities')
        rt=minmax_scale(rankdata(table.target_score,method='dense'))
        rg=minmax_scale(rankdata(table.G_score,method='dense'))
        final=minmax_scale(rankdata(rt+rg,method='dense'))
        for key,expected in [('target_score_scaled',rt),('G_score_scaled',rg),('influence_score',final)]:
            if not np.allclose(table[key],expected,rtol=1e-6,atol=1e-7):raise ValueError('Standard ANANSE component mismatch: '+key)
        de=pd.read_csv(contrast['path'],sep='\t').set_index('gene_id')
        if not (de.loc[table.factor,'log2FoldChange']>0).all():raise ValueError('Native positive-logFC eligibility violated')
        if sha(folder/'directed_de.tsv')!=sha(contrast['path']):raise ValueError('Directional DE checkpoint changed')
        reports.append(dict(mode=mode,source=source,target=target,ranking=str(ranking),diffnetwork=str(diffnet),de=contrast['path']))
        ranks.append(dict(source=source,target=target,TFs=len(table),top20=table.factor.head(20).tolist(),sha256=sha(ranking)))
    manifest=root/'report_qualification_manifest.tsv'
    with manifest.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,list(reports[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(reports)
    qualify(root,manifest)
    # Forward/reverse RNA scales and state means must be identical after swapping.
    if len(contrasts)!=2:raise ValueError('Expected both directions')
    heatmaps=[]
    for contrast in contrasts:
        plots=base/(contrast['source']+'_to_'+contrast['target'])/'main'/'plots'
        heatmaps.append(json.loads((plots/'heatmap_settings.json').read_text()))
    if any(heatmaps[0][key]!=heatmaps[1][key] for key in ['tpm_limits','fc_limits','width_mm','font_family']):raise ValueError('Contrast heatmap scales differ')
    receipt=dict(status='QUALIFIED_STATE_SPECIFIC_TPM_CACHED_AFFINITY',
        input_receipt=bound,expression_selection=expression,networks=networks,rankings=ranks,
        scope='Cached corrected-column-order affinities; standard multiplicative edges and native influence; no upstream recalculation or historical comparison',
        publication_status='NOT_REQUESTED' if not destination else 'PENDING')
    out=root/'qualification.json';out.write_text(json.dumps(receipt,indent=2)+'\n')
    if destination:
        dest=Path(destination)
        if dest.exists():raise ValueError('Refusing to overwrite an existing release')
        try:
            shutil.copytree(root,dest)
            # The receipt is finalized after copying; exclude its self-reference.
            copied={str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file() and p!=out}
            if any(sha(dest/relative)!=digest for relative,digest in copied.items()):raise ValueError('Published bytes differ')
            receipt['publication_status']='CHECKSUM_VERIFIED_COPY'
            receipt['published_destination']=str(dest)
            receipt['published_file_sha256']=copied
            out.write_text(json.dumps(receipt,indent=2)+'\n');shutil.copy2(out,dest/'qualification.json')
        except OSError as error:
            receipt['publication_status']='FAILED_PROJECT_COPY'
            receipt['publication_error']=str(error)
            out.write_text(json.dumps(receipt,indent=2)+'\n')
            raise
    print(json.dumps(dict(status=receipt['status'],publication=receipt['publication_status'])))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for key in ['results','params','input-receipt']:parser.add_argument('--'+key,required=True)
    parser.add_argument('--destination')
    run(**vars(parser.parse_args()))
