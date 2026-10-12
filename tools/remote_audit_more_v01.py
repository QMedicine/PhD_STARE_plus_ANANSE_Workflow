import pathlib,json,base64,hashlib
root=pathlib.Path('/vast/projects/Immune_Aging');paths=[]
tf=root/'.nextflow_assets/nf-core/tfactivity'
for folder in ['modules/local/counts/combine','modules/local/prepare_genome','subworkflows/local/prepare_genome']:
 p=tf/folder
 if p.is_dir():paths += [x for x in p.rglob('*') if x.is_file() and x.suffix in ['.py','.nf','.sh','.config']]
st=root/'Itokawa_et_al_2022_STARE/STARE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM'
paths += list(st.glob('*/*metadata*'))
an=root/'Itokawa_et_al_2022_ANANSE/STARE_output_in_ANANSE/STARE_output_in_ANANSE_results_100kb_w5000000_m5000000_H13_motifsV2_loess_offsets_CPM'
paths += list(an.glob('*/*influence.txt'))
paths += [root/'Itokawa_et_al_2022_ATACseq/sample_information.csv']
for p in sorted(set(paths)):
 if p.is_file() and p.stat().st_size<2000000:
  d=p.read_bytes();print(json.dumps({'path':str(p),'sha256':hashlib.sha256(d).hexdigest(),'content_base64':base64.b64encode(d).decode()}))
