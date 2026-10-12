import os,json,pathlib,hashlib,base64
root=pathlib.Path('/vast/projects/Immune_Aging')
needles=['004_HiC_analysis','005_HiC_analysis','006_HiC_analysis','007_HiC_analysis','01_Itokawa_et_al_2022_ATACseq_job_array.sh','02_Itokawa_et_al_2022_ATACseq_job_array_remove_blacklist_call_peaks.sh']
bases=[root/'20240603_MMP_Hi_Low_HSCs',root/'Itokawa_et_al_2022_ATACseq',root/'Itokawa_et_al_2022_nfcore_tfactivity',root/'Itokawa_et_al_2022_STARE']
for base in bases:
 print(json.dumps({'type':'directory','path':str(base),'exists':base.exists()}))
 if not base.exists(): continue
 for current,dirs,files in os.walk(base):
  depth=len(pathlib.Path(current).relative_to(base).parts)
  dirs[:]=sorted(d for d in dirs if d not in ['work','.git','.nextflow','fastq','fastqs','logs'] and depth<4)
  for f in sorted(files):
   p=pathlib.Path(current)/f
   wanted=any(n in f for n in needles) or (f.endswith(('.sh','.R','.py','.csv','.config','.yml','.yaml','.json')) and depth<2)
   if wanted and p.stat().st_size<500000:
    data=p.read_bytes(); print(json.dumps({'type':'script','path':str(p),'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'content_base64':base64.b64encode(data).decode()}))
   elif f.endswith(('.bam','.cool','.mcool','.h5','.narrowPeak','summits.bed','.transfac','.psem','.sif')):
    print(json.dumps({'type':'asset','path':str(p),'size':p.stat().st_size}))
