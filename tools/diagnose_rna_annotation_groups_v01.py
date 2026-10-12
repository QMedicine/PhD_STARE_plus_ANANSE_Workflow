#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd
project=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq');root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');out=root/'evidence/rna_discrepancy_v01'
old=pd.read_csv(project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv',index_col=0)
raw=pd.read_csv(project/'featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt',sep='\t',comment='#',index_col=0)
raw.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in raw.columns]
full=pd.read_csv(project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.tsv',sep='\t',index_col=0)
print('FULL_EXPORT_COLUMNS',list(full.columns))
summary={}
for style,names in [('raw_names',raw['gene_name']),('exported_names',full['gene_name'].reindex(raw.index))]:
 names=names.where(names.notna(),pd.Series(raw.index,index=raw.index));groups=names.groupby(names).groups
 grouped=raw[old.columns].groupby(names).sum()
 rebuilt=pd.DataFrame(grouped.loc[names.loc[old.index]].to_numpy(),index=old.index,columns=old.columns)
 mismatch=rebuilt.to_numpy()!=old.to_numpy();summary[style]=dict(n_different_values=int(mismatch.sum()))
 for i,j in zip(*np.where(mismatch)):
  gid=old.index[i];sample=old.columns[j];members=list(groups[names.loc[gid]])
  print('MISMATCH',style,gid,sample,'historical',old.iloc[i,j],'group_sum',rebuilt.iloc[i,j],'members',members,'values',list(raw.loc[members,sample]))
 if not mismatch.any():
  ledger=[dict(retained_gene_id=gid,member_gene_id=member,group_symbol=names.loc[gid],n_members=len(groups[names.loc[gid]]),style=style) for gid in old.index for member in groups[names.loc[gid]]]
  pd.DataFrame(ledger).to_csv(out/'authenticated_group_ledger.tsv',sep='\t',index=False)
(out/'annotation_group_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
for p in (project/'DEG_analysis_STAR_aligned_counts_all_cell_types').rglob('*'):
 if p.is_file() and p.suffix in ['.R','.Rmd','.py','.sh']:print('CLEANUP_CANDIDATE',p)
