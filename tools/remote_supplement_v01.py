import pathlib,json,base64,hashlib
root=pathlib.Path('/vast/projects/Immune_Aging'); paths=[]
hi=root/'20240603_MMP_Hi_Low_HSCs/HiCExplorer_analysis_of_OvsY_MMP_combined_HSCs/scripts'
paths+=list(hi.glob('*.sh'))
at=root/'Itokawa_et_al_2022_ATACseq'
paths += [at/'Itokawa_et_al_2022_ATACseq_job_array.sh',at/'Itokawa_et_al_2022_ATACseq_job_array_remove_blacklist_call_peaks.sh']
tf=root/'.nextflow_assets/nf-core/tfactivity'
for folder in ['subworkflows/local/counts','modules/local/counts/filter_genes','modules/local/counts/calculate_tpm','modules/local/prepare_genome','modules/nf-core/deseq2/differential','conf']:
 p=tf/folder
 if p.is_dir():paths+= [x for x in p.rglob('*') if x.is_file() and x.suffix in ['.nf','.py','.R','.config','.yml']]
paths += list((tf/'subworkflows/local').glob('*counts*/*.nf'))
for run in ['HSC_50kbp_window_H13CORE_motifs','HSC_50kbp_window_H13CORE_motifs_Young_over_Old']:
 p=root/'Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results'/run/'pipeline_info'
 paths+=list(p.glob('*.yml'))+list(p.glob('params*.json'))+list(p.glob('*trace*.txt'))
 paths+=list((p.parent/'00_prepare_genome').rglob('*length*.tsv'))
paths+=list((root/'Itokawa_et_al_2022_nfcore_tfactivity').glob('*HSC*.csv'))
an=pathlib.Path('/home/users/allstaff/stacpoole.q/.conda/envs/ananse/lib')
for p in an.glob('python*/site-packages/ananse'):
 paths += list(p.glob('*.py'))
for p in sorted(set(paths)):
 if p.is_file() and p.stat().st_size < 2000000:
  d=p.read_bytes();print(json.dumps({'path':str(p),'size':len(d),'sha256':hashlib.sha256(d).hexdigest(),'content_base64':base64.b64encode(d).decode()}))
