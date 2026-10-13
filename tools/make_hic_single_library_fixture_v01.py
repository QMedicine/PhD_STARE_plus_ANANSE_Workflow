#!/usr/bin/env python3
"""Create bounded paired prefixes of the declared practice inputs on scratch."""
import argparse,csv,gzip,hashlib,io,json,re
from pathlib import Path
import pysam

def build(samplesheet,outdir,records):
 out=Path(outdir)
 if out.exists():raise ValueError('Refuse to overwrite fixture inputs')
 out.mkdir(parents=True)
 with Path(samplesheet).open() as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
 audits=[]
 if 'pool_group' in fields:raise ValueError('This fixture must test a samplesheet without the optional pooling column')
 for row in rows:
  if row['assay']!='HiC':continue
  source=[Path(row['fastq_1']),Path(row['fastq_2'])];targets=[out/(row['sample_name']+f'_{mate}.fastq.gz') for mate in [1,2]]
  digest=hashlib.sha256()
  handles=[gzip.open(p,'wt') for p in targets]
  try:
   with pysam.FastxFile(str(source[0])) as f1,pysam.FastxFile(str(source[1])) as f2:
    for i in range(records):
     try:r1,r2=next(f1),next(f2)
     except StopIteration:raise ValueError('Practice library shorter than the requested prefix')
     if re.sub('/[12]$','',r1.name)!=re.sub('/[12]$','',r2.name):raise ValueError('Mate names differ')
     for handle,read in zip(handles,[r1,r2]):
      if read.quality is None or len(read.sequence)!=len(read.quality):raise ValueError('Invalid FASTQ quality')
      header=read.name+(' '+read.comment if read.comment else '')
      text='@'+header+'\n'+read.sequence+'\n+\n'+read.quality+'\n'
      handle.write(text);digest.update(text.encode())
  finally:
   for handle in handles:handle.close()
  audits.append(dict(sample_name=row['sample_name'],records=records,source_paths=list(map(str,source)),source_bytes=[p.stat().st_size for p in source],paired_prefix_sha256=digest.hexdigest(),fixture_paths=list(map(str,targets))))
  row['fastq_1'],row['fastq_2']=map(str,targets)
 with (out/'samplesheet.tsv').open('w') as f:
  writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
 (out/'fixture_provenance.json').write_text(json.dumps(dict(selection='first_N_paired_records_per_library',full_input_qualification=False,libraries=audits),indent=2)+'\n')
 return len(audits)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--samplesheet',required=True);p.add_argument('--outdir',required=True);p.add_argument('--records',type=int,default=5000)
 print(build(**vars(p.parse_args())))
