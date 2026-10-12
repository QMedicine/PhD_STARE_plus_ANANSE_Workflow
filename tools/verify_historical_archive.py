#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];checked=0;excluded=0
for manifest in ['historical/ananse_license_manifest_v01.json','historical/hic_gap_manifest_v01.json','historical/rna_cleanup_manifest_v01.json','historical/manifest_v01.json','historical/supplement_manifest_v01.json','historical/additional_manifest_v01.json','evidence/audit/discovery_manifest_v01.json']:
 for row in json.loads((root/manifest).read_text()):
  if 'local_path' not in row:continue
  path=root/row['local_path']
  if path.name.startswith('STARE_') and path.name.endswith('TF_influence.txt'):
   excluded+=1;continue
  if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise SystemExit('Checksum mismatch: '+str(path))
  checked+=1
print(json.dumps({'status':'PASS','verified_archive_copies':checked,'excluded_scientific_output_fixtures':excluded}))
