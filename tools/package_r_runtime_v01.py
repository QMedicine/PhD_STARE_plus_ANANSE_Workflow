#!/usr/bin/env python3
import csv,hashlib,json,os,re,subprocess,tarfile
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
with (root/'containers/r_library_closure_v01.tsv').open() as f:packages=list(csv.DictReader(f,delimiter='\t'))
archive=root/'containers/r_runtime_v01.tar'
if archive.exists():raise SystemExit('Refusing archive overwrite')
external={};record=[]
with tarfile.open(archive,'w') as tar:
    for row in packages:
        p=Path(row['LibPath'])/row['Package']
        tar.add(p,arcname='Rlibrary/'+row['Package'],recursive=True)
        for so in p.rglob('*.so'):
            result=subprocess.run(['ldd',str(so)],text=True,capture_output=True,check=True)
            for line in result.stdout.splitlines():
                if 'not found' in line:raise ValueError('Unresolved native library '+str(so)+' '+line)
                m=re.search(r'=> (/\S+)',line)
                if not m:continue
                lib=Path(m[1]);name=lib.name
                if name.startswith(('libc.so','libm.so','libpthread.so','librt.so','libdl.so','libR.so','ld-linux')):continue
                digest=hashlib.sha256(lib.read_bytes()).hexdigest()
                if name in external and external[name]['sha256']!=digest:raise ValueError('Conflicting library basename '+name)
                external[name]=dict(path=str(lib),sha256=digest)
    for name,r in sorted(external.items()):tar.add(Path(r['path']).resolve(),arcname='host-libs/'+name)
(root/'evidence/r_runtime_v01_dependencies.json').write_text(json.dumps(dict(packages=packages,native_libraries=external),indent=2)+'\n')
