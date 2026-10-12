#!/usr/bin/env python3
import csv,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
project='/vast/projects/Immune_Aging/';scratch='/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow';run=project+'PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs'
tf=project+'Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/'
assets=list(map(json.loads,(root/'evidence/audit/curated_assets_observed_v01.jsonl').read_text().splitlines()))
roles={'fasta':'fasta','gtf':'gtf','chrom_sizes':'chrom_sizes','blacklist_atac':'atac_blacklist','blacklist_summits':'summit_blacklist','cleaned_counts':'selected_counts','tpm':'tpm','idmap':'idmap','psem':'psem','stare_sif':'stare_container','deseq_sif':'deseq2_container'}
common={roles[r['role']]:r['path'] for r in assets if r['role'] in roles}
common.update(atac_modes=['historical_column_order','corrected_column_order'],historical_header_samples=[r['sample'] for age in ['Young','Old'] for r in assets if r['role']=='bam' and r['stage']=='ATAC' and '_'+age+'_' in r['sample']],historical_data_samples=[r['sample'] for age in ['Old','Young'] for r in assets if r['role']=='bam' and r['stage']=='ATAC' and '_'+age+'_' in r['sample']],deseq_state_order=['Old_HSC','Young_HSC'],container_binds=['/vast/projects/Immune_Aging',scratch],samplesheet=run+'/samplesheet_Old_and_Young_HSCs_v01.tsv',source2target=run+'/source2target_Old_and_Young_HSCs_v01.tsv',reverse_contrast=True,gene_lengths=tf+'HSC_50kbp_window_H13CORE_motifs/00_prepare_genome/04_gtftools_length/gtf.txt',ananse_container=scratch+'/containers/ananse_v03.sif',r_plot_container=scratch+'/containers/r_plot_v01.sif',hicexplorer_container=scratch+'/containers/hicexplorer_v01.sif',tpm_container=project+'apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img',featurecounts_container=scratch+'/containers/bam_tools_v01.sif',atac_tools_container=scratch+'/containers/bam_tools_v01.sif',bedtools_container=project+'apptainer_singularity_cache/depot.galaxyproject.org-singularity-bedtools-2.31.1--hf5e1c6e_0.img')
for filename,role in [('affinities','historical_affinity'),('hic_checkpoints','historical_corrected_h5')]:
 with (root/f'templates/{filename}_Old_and_Young_HSCs_v01.tsv').open('w') as f:
  w=csv.DictWriter(f,['cell_state','path'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(dict(cell_state=r['sample'],path=r['path']) for r in assets if r['role']==role)
rows=[dict(source='Old_HSC',target='Young_HSC',path=tf+'HSC_50kbp_window_H13CORE_motifs_Young_over_Old/01_counts/05_deseq2/A_Old_HSC:B_Young_HSC/A_Old_HSC:B_Young_HSC.deseq2.results.tsv'),dict(source='Young_HSC',target='Old_HSC',path=tf+'HSC_50kbp_window_H13CORE_motifs/01_counts/05_deseq2/A_Young_HSC:B_Old_HSC/A_Young_HSC:B_Old_HSC.deseq2.results.tsv')]
with (root/'templates/de_Old_and_Young_HSCs_v01.tsv').open('w') as f:
 w=csv.DictWriter(f,['source','target','path'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
for stage in ['preflight','rna_recount','historical_replay','bam_historical_hic_checkpoint','atac_historical_hic_checkpoint']:
 p=dict(common,stage=stage)
 if stage=='historical_replay':p.update(affinities=run+'/affinities_Old_and_Young_HSCs_v01.tsv',de_manifest=run+'/de_Old_and_Young_HSCs_v01.tsv')
 if stage in ['bam_historical_hic_checkpoint','atac_historical_hic_checkpoint']:p.update(hic_checkpoints=run+'/hic_checkpoints_Old_and_Young_HSCs_v01.tsv')
 if stage=='atac_historical_hic_checkpoint':p.update(de_manifest=run+'/de_Old_and_Young_HSCs_v01.tsv')
 (root/f'templates/params_{stage}_v01.json').write_text(json.dumps(p,indent=2)+'\n')
