#!/usr/bin/env python3
import hashlib,json,re,shutil,subprocess,tarfile
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
p=Path(shutil.which('ps')).resolve();archive=root/'containers/procps_v01.tar'
if archive.exists():raise SystemExit('Refusing archive overwrite')
libs={}
result=subprocess.run(['ldd',str(p)],capture_output=True,text=True,check=True)
if 'not found' in result.stdout:raise ValueError('ps dependency missing')
for line in result.stdout.splitlines():
 m=re.search(r'=> (/\S+)',line)
 if not m:continue
 lib=Path(m[1]);leaf=lib.name
 if leaf.startswith(('libc.so','libm.so','libpthread.so','librt.so','libdl.so','ld-linux')):continue
 libs[leaf]=dict(source=str(lib),sha256=hashlib.sha256(lib.read_bytes()).hexdigest())
with tarfile.open(archive,'w') as tar:
 tar.add(p,arcname='bin/ps')
 for leaf,entry in libs.items():tar.add(Path(entry['source']).resolve(),arcname='host-libs/'+leaf)
(root/'evidence/procps_v01_inventory.json').write_text(json.dumps(dict(binary_source=str(p),binary_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),libraries=libs),indent=2)+'\n')
