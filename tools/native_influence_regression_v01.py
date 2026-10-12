#!/usr/bin/env python3
"""Compare authentic saved influence components and edges with native replay outputs."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
original=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM')
records=[]
for source,target,direction,stem in [('Old_HSC','Young_HSC','Old2Young','STARE_OldHSC_to_YoungHSC_TF_influence'),('Young_HSC','Old_HSC','Young2Old','STARE_YoungHSC_to_OldHSC_TF_influence')]:
 expected=original/('TF_Influence_'+direction)/(stem+'.txt')
 observed=root/'runs/network_replay_v09/results/archived_affinity_replay'/(source+'_to_'+target)/'tables/influence.txt'
 a=pd.read_csv(expected,sep='\t',index_col='factor');b=pd.read_csv(observed,sep='\t',index_col='factor')
 if set(a.index)!=set(b.index) or list(a.columns)!=list(b.columns):raise ValueError('Influence identities/component columns differ')
 aligned=b.loc[a.index,a.columns];equal=(a.to_numpy()==aligned.to_numpy())|(pd.isna(a).to_numpy()&pd.isna(aligned).to_numpy())
 delta=np.abs(a.to_numpy(dtype=float)-aligned.to_numpy(dtype=float))
 ea=original/('TF_Influence_'+direction)/(stem+'_diffnetwork.tsv');eb=observed.parent/'influence_diffnetwork.tsv'
 x=pd.read_csv(ea,sep='\t').set_index(['source','target']);y=pd.read_csv(eb,sep='\t').set_index(['source','target'])
 if x.index.has_duplicates or y.index.has_duplicates:raise ValueError('Duplicate differential-network edges')
 same_edges=set(x.index)==set(y.index);weights_equal=False
 if same_edges:weights_equal=bool((x['weight'].to_numpy()==y.loc[x.index,'weight'].to_numpy()).all())
 record=dict(source=source,target=target,ranking_same_factor_order=list(a.index)==list(b.index),ranking_byte_match=expected.read_bytes()==observed.read_bytes(),ranking_components_exact=bool(equal.all()),different_component_values=int((~equal).sum()),max_absolute_component_error=float(np.nanmax(delta)),top20_same_factors=set(a.head(20).index)==set(b.head(20).index),differential_same_edges=same_edges,differential_weights_exact=weights_equal,differential_byte_match=ea.read_bytes()==eb.read_bytes(),expected_sha256=hashlib.sha256(expected.read_bytes()).hexdigest(),observed_sha256=hashlib.sha256(observed.read_bytes()).hexdigest())
 records.append(record)
Path('influence_regression.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
if not all(r['ranking_components_exact'] and r['differential_same_edges'] and r['differential_weights_exact'] for r in records):raise SystemExit('Influence regression differs; do not declare exact reproduction')
