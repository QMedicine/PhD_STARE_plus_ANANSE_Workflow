#!/usr/bin/env python3
"""Validate actual STARE products and biological-pool count conservation in fixtures."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
folder=root/'runs/fixed_contacts_stare_smoke_v28/results/diffhic_loess_fixed_grid__synthetic_smoke/stare'
records=[]
for state in ['State_1','State_2']:
    files=list((folder/state/state/'Gene_TF_matrices').glob('*TF_Gene_Affinities*.txt.gz'))
    if len(files)!=2:raise ValueError('Expected two replicate STARE affinity matrices: '+str(files))
    for path in files:
        matrix=pd.read_csv(path,sep='\t',index_col=0)
        values=matrix.to_numpy(dtype=float)
        if set(matrix.index)!=set(['SYNTH_GENE_1','SYNTH_GENE_2','SYNTH_GENE_3']) or matrix.shape[1]!=759 or not np.isfinite(values).all() or not (values>0).any():raise ValueError('Synthetic STARE matrix identities or scores differ: '+str((path,matrix.shape,matrix.index.tolist())))
        records.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),genes=len(matrix),numeric_columns=matrix.shape[1],finite=True,positive_values=int((values>0).sum())))
pool=root/'runs/hic_pooling_fixture_v27/results/diffhic/fixed_grid_contacts'
counts=pd.read_csv(pool/'fixed_grid_counts.tsv',sep='\t')
labels=['State_1_POOL','State_2_POOL']
if counts[labels].sum().tolist()!=[14,14]:raise ValueError('Pooling changed biological library count conservation')
provenance=json.loads((pool/'fixed_contacts_provenance.json').read_text())
if provenance['library_order']!=labels or provenance['fixed_grid_count_totals']!=[14,14]:raise ValueError('Pooled identities/totals differ')
out=root/'evidence/connector_fixture_v01';out.mkdir(exist_ok=True)
receipt=dict(status='PASS_FIXED_CONTACTS_TO_STARE_AND_POOLING',stare_affinity_matrices=records,pooling=dict(components=4,biological_pools=2,retained_pairs_per_pool=[14,14]),qualification_boundary='Synthetic native integration only; full-cohort raw Hi-C and original custom-mapper equivalence not claimed')
(out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
