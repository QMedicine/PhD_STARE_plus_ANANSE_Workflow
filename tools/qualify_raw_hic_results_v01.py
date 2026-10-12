#!/usr/bin/env python3
"""Qualify and checksum-copy new-method GRNs/tables without claiming historical identity."""
import csv,hashlib,json,shutil
from pathlib import Path
import numpy as np,pandas as pd
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');run=root/'runs/hic_raw_to_stare_v29/results'
project=Path('/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/validated_results/hic_raw_to_stare_v29')
qual=list((root/'evidence/raw_hic_bam_qualification_v01').glob('*/qualification.json'))
if len(qual)!=12 or any(json.loads(p.read_text())['status']!='QUALIFIED_RAW_HIC_BAM_INTERNAL_CONTRACT' for p in qual):raise ValueError('All 12 raw BAMs must qualify before result publication')
contact=json.loads((run/'diffhic/fixed_grid_contacts/fixed_contacts_provenance.json').read_text())
if len(contact['library_order'])!=6 or len(contact['states'])!=2 or contact['method']!='diffhic_loess_fixed_grid':raise ValueError('Full-cohort pooling/contact identities differ')
scores=[];networks=[]
for mode in ['historical_column_order','corrected_column_order']:
    path=run/('diffhic_loess_fixed_grid__'+mode)
    files=list(path.glob('networks/*/*_SYMBOLS.tsv'))
    for network in files:
        n=0
        for table in pd.read_csv(network,sep='\t',chunksize=500000):
            if table.columns.tolist()!=['tf_target','prob'] or not np.isfinite(table['prob']).all() or not table['prob'].between(0,1).all() or table['tf_target'].isna().any():raise ValueError('Invalid network schema/probability')
            n+=len(table)
        if n==0:raise ValueError('Empty GRN')
        networks.append(dict(mode=mode,path=str(network.relative_to(run)),edges=n))
    ranking_files=list(path.glob('*_to_*/tables/influence.txt'))
    if len(ranking_files)!=2:raise ValueError('Expected two directed influence products per ATAC mode')
    for ranking in ranking_files:
        table=pd.read_csv(ranking,sep='\t');components=['influence_score','target_score_scaled','G_score_scaled']
        if len(table)<20 or table['factor'].duplicated().any() or not np.isfinite(table[components]).all().all() or not table[components].apply(lambda col:col.between(0,1).all()).all():raise ValueError('Influence identity/scaling validation failed')
        scores.append(dict(mode=mode,path=str(ranking.relative_to(run)),tfs=len(table),top20=table['factor'].head(20).tolist()))
if len(networks)!=4:raise ValueError('Expected four state GRNs across two ATAC modes')
allowed=[]
for p in run.rglob('*'):
    if not p.is_file():continue
    relative=p.relative_to(run)
    final_mode_product=('diffhic_loess_fixed_grid__' in str(relative) and ('networks' in relative.parts or 'tables' in relative.parts))
    count_product=('diffhic' in relative.parts and p.suffix in ['.tsv','.json','.txt'])
    if p.suffix in ['.tsv','.txt','.json'] and (final_mode_product or 'provenance' in relative.parts or count_product):
        allowed.append((p,relative))
project.mkdir(parents=True,exist_ok=True);copies=[]
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
for p,relative in allowed:
    target=project/relative;target.parent.mkdir(parents=True,exist_ok=True)
    sha=digest(p)
    if target.exists():
        if digest(target)!=sha:raise ValueError('Refuse differing durable result overwrite')
    else:shutil.copy2(p,target)
    if digest(target)!=sha:raise ValueError('Durable result copy differs')
    copies.append(dict(path=str(relative),bytes=p.stat().st_size,sha256=sha))
receipt=dict(status='QUALIFIED_RAW_HIC_INTERNAL_CONTRACT_GRNS_AND_TABLES_REPORTS_INCOMPLETE',method='diffhic_loess_fixed_grid',raw_libraries=12,biological_pools=6,states=2,networks=networks,influence=scores,verified_copies=copies,historical_t32_equivalence='UNKNOWN',historical_ICE_equivalence='NOT_CLAIMED',reports='WITHHELD_PENDING_USER_PLOT_DECISION')
(project/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n')
out=root/'evidence/raw_hic_results_v01';out.mkdir(exist_ok=True);(out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='verified_copies'},indent=2))
