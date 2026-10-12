import pathlib,os,json,base64,hashlib
root=pathlib.Path('/vast/projects/Immune_Aging')
for base in [root/'Itokawa_et_al_2022_RNAseq',root/'Itokawa_et_al_2022_ATACseq',root/'20240603_MMP_Hi_Low_HSCs/raw_fastq_MMP_high_low',root/"20240603_MMP_Hi_Low_HSCs/Hannah's_Hi-C_fastq_alignment_QC_processing_scripts"]:
 if not base.exists():print(json.dumps(dict(type='directory',path=str(base),exists=False)));continue
 for current,dirs,files in os.walk(base):
  depth=len(pathlib.Path(current).relative_to(base).parts)
  dirs[:]=sorted(d for d in dirs if depth<2 and d not in ['STAR_aligned_bams','featureCounts','filtered','macs3_peaks','DEG_analysis_STAR_aligned_counts_all_cell_types','qc_reports','logs','logs_QC'])
  for name in sorted(files):
   p=pathlib.Path(current)/name
   if name.endswith(('.fastq.gz','.fq.gz','.fastq','.fq')):
    print(json.dumps(dict(type='fastq',path=str(p),size=p.stat().st_size)))
   elif p.suffix in ['.sh','.slurm','.sbatch','.config','.csv','.tsv'] and p.stat().st_size<200000 and (depth==0 or 'Hannah' in str(base)):
    d=p.read_bytes();print(json.dumps(dict(type='supporting_script',path=str(p),size=len(d),sha256=hashlib.sha256(d).hexdigest(),content_base64=base64.b64encode(d).decode())))
