import json
from pathlib import Path
import pandas as pd,numpy as np
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
original=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM')
records=[]
for source,target,direction,stem in [('Old_HSC','Young_HSC','Old2Young','STARE_OldHSC_to_YoungHSC_TF_influence'),('Young_HSC','Old_HSC','Young2Old','STARE_YoungHSC_to_OldHSC_TF_influence')]:
 a=pd.read_csv(original/('TF_Influence_'+direction)/(stem+'.txt'),sep='\t',index_col='factor')
 b=pd.read_csv(root/'runs/network_replay_v09/results/archived_affinity_replay'/(source+'_to_'+target)/'tables/influence.txt',sep='\t',index_col='factor')
 cols={c:int((a[c]!=b.loc[a.index,c]).sum()) for c in a.columns}
 identical_score_sequence=np.array_equal(a.influence_score.to_numpy(),b.influence_score.to_numpy())
 changes=[dict(position=i+1,expected=str(tf),observed=str(b.index[i]),expected_score=float(a.iloc[i].influence_score),observed_score=float(b.iloc[i].influence_score)) for i,tf in enumerate(a.index) if b.index[i]!=tf]
 records.append(dict(source=source,target=target,component_differences=cols,top20_order_identical=list(a.head(20).index)==list(b.head(20).index),same_score_sequence=identical_score_sequence,changed_positions=changes))
print(json.dumps(records,indent=2));(root/'evidence/influence_regression_v01/score_order_detail.json').write_text(json.dumps(records,indent=2)+'\n')
