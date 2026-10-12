#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import pandas as pd
project=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq')
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
selection=project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv'
old=pd.read_csv(selection,index_col=0)
new=pd.read_csv(root/'runs/rna_recount_v01/results/rna/recount/selected_counts.csv',index_col=0)
raw=pd.read_csv(project/'featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt',sep='\t',comment='#',index_col=0)
raw.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in raw.columns]
new=new.loc[old.index,old.columns];raw_values=raw.loc[old.index,old.columns]
rows=[]
for gid,sample in zip(*np.where(old.to_numpy()!=new.to_numpy())):
 gene=old.index[gid];column=old.columns[sample]
 rows.append(dict(gene_id=gene,sample=column,historical_cleaned=old.iloc[gid,sample],original_raw_featurecounts=raw_values.loc[gene,column],recount=new.iloc[gid,sample],gene_name=raw.loc[gene,'gene_name'],gene_biotype=raw.loc[gene,'gene_biotype']))
out=root/'evidence/rna_discrepancy_v01';out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out/'mismatches.tsv',sep='\t',index=False)
summary=dict(n_genes=len(old),n_samples=len(old.columns),cleaned_vs_recount_differences=len(rows),raw_vs_recount_differences=int((raw_values.to_numpy()!=new.to_numpy()).sum()),raw_vs_cleaned_differences=int((raw_values.to_numpy()!=old.to_numpy()).sum()),affected_genes=sorted({r['gene_id'] for r in rows}))
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
print(pd.DataFrame(rows).head(15).to_string(index=False))
