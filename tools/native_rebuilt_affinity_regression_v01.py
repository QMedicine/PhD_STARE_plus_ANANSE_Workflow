#!/usr/bin/env python3
"""Test each explicitly labelled ATAC ordering mode against saved affinities and GRNs."""
import csv,gzip,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
original=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/STARE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM')
original_network=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM')
records=[];digest_cache={}
for mode in ['historical_column_order','corrected_column_order']:
 for state,old_label,age in [('Old_HSC','Old_LTHSC','Old'),('Young_HSC','Young_LTHSC','Young')]:
  observed_dir=root/'runs/atac_checkpoint_v05/results'/mode/'stare'/state
  for replicate in [1,2]:
   expected=original/old_label/'Gene_TF_matrices'/f'{old_label}_TF_Gene_Affinities_{age}_rep{replicate}.txt.gz'
   files=list(observed_dir.rglob(f'*TF_Gene_Affinities_*Rep{replicate}.txt.gz'))
   if len(files)!=1:raise ValueError('Expected exactly one reconstructed affinity file')
   observed=files[0]
   a=pd.read_csv(expected,sep='\t',index_col=0);b=pd.read_csv(observed,sep='\t',index_col=0)
   same_keys=set(a.index)==set(b.index);same_columns=set(a.columns)==set(b.columns)
   record=dict(mode=mode,cell_state=state,replicate=replicate,expected_shape=list(a.shape),observed_shape=list(b.shape),same_gene_ids=same_keys,same_columns=same_columns)
   if same_keys and same_columns:
    b=b.loc[a.index,a.columns]
    numeric=[c for c in a.columns if pd.api.types.is_numeric_dtype(a[c])]
    av=a[numeric].to_numpy(dtype=float);bv=b[numeric].to_numpy(dtype=float)
    equal=(av==bv)|(np.isnan(av)&np.isnan(bv));finite=np.isfinite(av)&np.isfinite(bv)
    record.update(numeric_columns=len(numeric),different_values=int((~equal).sum()),numeric_exact=bool(equal.all()),max_absolute_error=float(np.abs(av[finite]-bv[finite]).max(initial=0)),within_1e_12=bool(np.allclose(av,bv,rtol=1e-12,atol=1e-12,equal_nan=True)))
   records.append(record);print(json.dumps(record),flush=True)
  expected_network=original_network/(state+'_STARE_to_ananse_network_SYMBOLS.tsv')
  observed_network=root/'runs/atac_checkpoint_v05/results'/mode/'networks'/state/(state+'_SYMBOLS.tsv')
  def sha(p):
   if str(p) in digest_cache:return digest_cache[str(p)]
   h=hashlib.sha256()
   with p.open('rb') as f:
    for block in iter(lambda:f.read(1048576),b''):h.update(block)
   digest_cache[str(p)]=h.hexdigest();return digest_cache[str(p)]
  records.append(dict(mode=mode,cell_state=state,network_byte_match=sha(expected_network)==sha(observed_network),expected_network_sha256=sha(expected_network),observed_network_sha256=sha(observed_network)))
Path('rebuilt_affinity_regression.json').write_text(json.dumps(records,indent=2)+'\n')
