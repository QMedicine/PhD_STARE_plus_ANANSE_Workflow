#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd
project=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq');root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
old=pd.read_csv(project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv',index_col=0)
raw=pd.read_csv(project/'featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt',sep='\t',comment='#',index_col=0)
raw.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in raw.columns]
values=raw[old.columns]
summary={};ledger=[]
for style in ['exact_symbol','uppercase_symbol']:
 names=raw['gene_name'].copy()
 names=names.where(names.notna(),pd.Series(raw.index,index=raw.index))
 if style=='uppercase_symbol':names=names.str.upper()
 grouped=values.groupby(names).sum()
 reconstructed=pd.DataFrame(grouped.loc[names.loc[old.index]].to_numpy(),index=old.index,columns=old.columns)
 mismatch=int((reconstructed.to_numpy()!=old.to_numpy()).sum());summary[style]=dict(n_different_values=mismatch,exact_numeric=mismatch==0)
 if mismatch==0:
  for gid in old.index:
   members=list(raw.index[names==names.loc[gid]])
   for member in members:ledger.append(dict(retained_gene_id=gid,member_gene_id=member,group_symbol=names.loc[gid],n_members=len(members),style=style))
  pd.DataFrame(ledger).to_csv(root/'evidence/rna_discrepancy_v01/duplicate_group_ledger.tsv',sep='\t',index=False)
  break
(root/'evidence/rna_discrepancy_v01/duplicate_group_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
for p in (project/'DEG_analysis_STAR_aligned_counts_all_cell_types').rglob('*'):
 if p.is_file() and p.suffix in ['.R','.Rmd','.py','.sh'] and ('dup' in p.name.lower() or 'clean' in p.name.lower()):print('CLEANUP_SCRIPT',p)
