#!/usr/bin/env python3
import hashlib,json,re,shutil,subprocess,tarfile
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');archive=root/'containers/bam_tools_v01.tar'
if archive.exists():raise SystemExit('Refusing archive overwrite')
libs={};record={}
with tarfile.open(archive,'w') as tar:
 for name in ['samtools','bedtools','featureCounts']:
  binary=Path(shutil.which(name)).resolve();record[name]=dict(source=str(binary),sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
  tar.add(binary,arcname='tools/'+name)
  result=subprocess.run(['ldd',str(binary)],capture_output=True,text=True,check=True)
  for line in result.stdout.splitlines():
   if 'not found' in line:raise ValueError('Unresolved dependency '+line)
   m=re.search(r'=> (/\S+)',line)
   if not m:continue
   lib=Path(m[1]);leaf=lib.name
   if leaf.startswith(('libc.so','libm.so','libpthread.so','librt.so','libdl.so','ld-linux')):continue
   sha=hashlib.sha256(lib.read_bytes()).hexdigest()
   if leaf in libs and libs[leaf]['sha256']!=sha:raise ValueError('Conflicting dependency '+leaf)
   libs[leaf]=dict(source=str(lib),sha256=sha)
 for name,data in libs.items():tar.add(Path(data['source']).resolve(),arcname='tool-libs/'+name)
(root/'evidence/bam_tools_v01_inventory.json').write_text(json.dumps(dict(binaries=record,libraries=libs),indent=2)+'\n')
