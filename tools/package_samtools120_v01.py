#!/usr/bin/env python3
import hashlib,json,re,shutil,subprocess,tarfile
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');archive=root/'containers/samtools120_v01.tar'
if archive.exists():raise SystemExit('Refuse archive overwrite')
binary=Path(shutil.which('samtools')).resolve()
version=subprocess.check_output([str(binary),'--version'],text=True)
if not version.startswith('samtools 1.20\nUsing htslib 1.20\n'):raise ValueError('Expected authenticated original samtools/htslib 1.20')
result=subprocess.run(['ldd',str(binary)],capture_output=True,text=True,check=True);libraries={}
for line in result.stdout.splitlines():
 if 'not found' in line:raise ValueError('Unresolved dependency '+line)
 match=re.search(r'=> (/\S+)',line)
 if not match:continue
 path=Path(match[1])
 if path.name.startswith(('libc.so','libm.so','libpthread.so','librt.so','libdl.so','ld-linux')):continue
 libraries[path.name]=dict(source=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
with tarfile.open(archive,'w') as tar:
 tar.add(binary,arcname='samtools')
 for name,data in libraries.items():tar.add(Path(data['source']),arcname='lib/'+name)
(root/'evidence/samtools120_v01_inventory.json').write_text(json.dumps(dict(binary=str(binary),sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),version=version,libraries=libraries),indent=2)+'\n')
