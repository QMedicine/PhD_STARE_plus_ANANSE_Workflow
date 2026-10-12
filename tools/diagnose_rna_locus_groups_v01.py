#!/usr/bin/env python3
"""Qualify the newly recovered same-chromosome/strand duplicate aggregation rule."""
import json
from pathlib import Path
import pandas as pd
project=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq')
raw=pd.read_csv(project/'featureCounts/gene_counts_HSCs_MPPs_etc_with_gene_names_biotype.txt',sep='\t',comment='#',index_col=0)
raw.columns=[Path(c).name.removesuffix('_Aligned.sortedByCoord.out.bam') for c in raw.columns]
selected=pd.read_csv(project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Counts_tables_no_ENSMUS_in_gene_name/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final_gene_name_dups_cleaned_no_ENSMUS_in_gene_name_HSC_Geneid.csv',index_col=0)
full=pd.read_csv(project/'DEG_analysis_STAR_aligned_counts_all_cell_types/Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.tsv',sep='\t',index_col=0)
groups=raw.groupby(full['gene_name'].reindex(raw.index),sort=False).groups
ledger=[];rebuilt={}
for gid in selected.index:
 name=full.loc[gid,'gene_name'];members=list(groups[name])
 chromosomes={x for value in raw.loc[members,'Chr'].astype(str) for x in value.split(';') if x}
 strands={x for value in raw.loc[members,'Strand'].astype(str) for x in value.split(';') if x}
 same_locus=len(chromosomes)==1 and len(strands)==1
 included=members if same_locus else [gid]
 rebuilt[gid]=raw.loc[included,selected.columns].sum()
 for member in included:ledger.append(dict(retained_gene_id=gid,member_gene_id=member,group_symbol=name,action='sum_same_chromosome_and_strand' if same_locus and len(members)>1 else 'retain_selected_representative',n_members=len(included)))
values=pd.DataFrame.from_dict(rebuilt,orient='index').loc[selected.index,selected.columns]
diff=values.to_numpy()!=selected.to_numpy()
summary=dict(n_genes=len(values),n_samples=len(selected.columns),n_different_values=int(diff.sum()),exact_numeric_match=not bool(diff.any()),summed_groups=len({r['retained_gene_id'] for r in ledger if r['n_members']>1}),selection_policy='Historical selected IDs are an explicit checkpoint; same chromosome/strand groups summed; distinct loci retain selected ID')
Path('locus_groups_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if not diff.any():pd.DataFrame(ledger).to_csv('historical_cleanup_ledger.tsv',sep='\t',index=False)
else:
 import numpy as np
 for i,j in zip(*np.where(diff)):print('MISMATCH',values.index[i],values.columns[j],values.iloc[i,j],selected.iloc[i,j])
