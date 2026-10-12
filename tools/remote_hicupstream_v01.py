import pathlib,json,base64,hashlib
p=pathlib.Path("/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/Hannah's_Hi-C_fastq_alignment_QC_processing_scripts/20240620_Hannah's_Processing_data_quality_control.Rmd")
d=p.read_bytes();print(json.dumps(dict(path=str(p),sha256=hashlib.sha256(d).hexdigest(),content_base64=base64.b64encode(d).decode())))
for base in [pathlib.Path('/vast/projects/Immune_Aging/.nextflow_assets/nf-core/tfactivity'),pathlib.Path('/home/users/allstaff/stacpoole.q/.conda/envs/ananse/lib/python3.10/site-packages/ananse-0.5.1.dist-info')]:
 for p in base.glob('*LICENSE*'):
  if p.is_file():
   d=p.read_bytes();print(json.dumps(dict(path=str(p),sha256=hashlib.sha256(d).hexdigest(),content_base64=base64.b64encode(d).decode())))
