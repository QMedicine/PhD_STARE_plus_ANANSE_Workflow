#!/usr/bin/env python3
"""Qualify and publish the independent-library practice GRNs and score tables."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path
import numpy as np,pandas as pd

def sha(path):
 digest=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):digest.update(chunk)
 return digest.hexdigest()
def qualify(results,samplesheet,bam_qualifications,destination):
 root=Path(results);dest=Path(destination)
 with Path(samplesheet).open() as f:reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;hic=[r for r in reader if r['assay']=='HiC']
 if 'pool_group' in fields:raise ValueError('This qualification expects a samplesheet without a pooling column')
 labels=[r['sample_name'] for r in hic]
 for row in hic:
  record=json.loads((Path(bam_qualifications)/row['sample_name']/'qualification.json').read_text())
  if record['status']!='QUALIFIED_RAW_HIC_BAM_INTERNAL_CONTRACT' or record['sample']!=row:raise ValueError('Original-input and qualified-BAM identities differ')
 contact=json.loads((root/'diffhic/fixed_grid_contacts/fixed_contacts_provenance.json').read_text())
 if contact['library_order']!=labels or set(contact['states'])!={r['cell_state'] for r in hic}:raise ValueError('Counting changed library/state identities')
 manifest=json.loads((root/'provenance/validated_inputs.json').read_text())
 if any(g['merge_required'] or g['routing']!='single_library' for g in manifest['hic_library_groups']):raise ValueError('Unexpected pooling route')
 for trace in (root/'provenance').glob('trace_*.tsv'):
  with trace.open() as f:
   if any('DIFFHIC_POOL_LIBRARIES' in r['name'] for r in csv.DictReader(f,delimiter='\t')):raise ValueError('A library merge task ran')
 networks=[];rankings=[]
 modes=['diffhic_loess_fixed_grid__historical_column_order','diffhic_loess_fixed_grid__corrected_column_order']
 for mode in modes:
  files=list((root/mode).glob('networks/*/*_SYMBOLS.tsv'))
  if len(files)!=len(contact['states']):raise ValueError('Wrong number of state GRNs')
  for file in files:
   edges=0
   for table in pd.read_csv(file,sep='\t',chunksize=500000):
    if list(table.columns)!=['tf_target','prob'] or table.tf_target.isna().any() or not np.isfinite(table.prob).all() or not table.prob.between(0,1).all():raise ValueError('Invalid GRN probabilities')
    edges+=len(table)
   if edges==0:raise ValueError('Empty GRN')
   networks.append(dict(mode=mode,path=str(file.relative_to(root)),edges=edges))
  files=list((root/mode).glob('*_to_*/tables/influence.txt'))
  if len(files)!=2:raise ValueError('Expected both influence directions')
  for file in files:
   table=pd.read_csv(file,sep='\t');scores=table[['influence_score','target_score_scaled','G_score_scaled']]
   if len(table)<20 or table.factor.duplicated().any() or not np.isfinite(scores).all().all() or not scores.apply(lambda x:x.between(0,1).all()).all():raise ValueError('Invalid influence table')
   rankings.append(dict(mode=mode,path=str(file.relative_to(root)),tfs=len(table)))
 if dest.exists():raise ValueError('Refuse to overwrite an existing qualified release')
 copies=[]
 for p in root.rglob('*'):
  if not p.is_file() or p.suffix not in ['.tsv','.txt','.json']:continue
  rel=p.relative_to(root)
  is_final=rel.parts[0] in modes and ('networks' in rel.parts or 'tables' in rel.parts)
  if not (is_final or 'provenance' in rel.parts or 'diffhic' in rel.parts):continue
  target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
  digest=sha(p)
  if digest!=sha(target):raise ValueError('Publication checksum differs')
  copies.append(dict(path=str(rel),sha256=digest,bytes=p.stat().st_size))
 receipt=dict(status='QUALIFIED_SINGLE_LIBRARY_HIC_GRNS_TABLES_REPORTS_INCOMPLETE',method='diffhic_loess_fixed_grid',raw_libraries=len(labels),library_names=labels,merge_tasks=0,duplicate_marking_scope='each supplied complete library',networks=networks,influence=rankings,reports='QUEUED_REPORT_ONLY_PHASE',verified_copies=copies)
 (root/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');shutil.copy2(root/'qualification.json',dest/'qualification.json');return receipt
if __name__=='__main__':
 p=argparse.ArgumentParser()
 for key in ['results','samplesheet','bam-qualifications','destination']:p.add_argument('--'+key,required=True)
 print(json.dumps(qualify(**vars(p.parse_args())),indent=2))
