#!/usr/bin/env python3
"""Snapshot observed module payloads; never alter their installed prefixes."""
import hashlib,json,re,shutil,subprocess,tarfile,os
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
archive=root/'containers/raw_tools_v02.tar'
if archive.exists():raise SystemExit('Refusing archive overwrite')
libs={};inventory=[];seen=set();vendor_dirs=[]
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def native_libraries(p):
 if p.is_symlink() or not p.is_file():return
 with p.open('rb') as f:magic=f.read(4)
 if magic!=b'\x7fELF':return
 env=dict(os.environ);env['LD_LIBRARY_PATH']=':'.join(vendor_dirs+[env.get('LD_LIBRARY_PATH','')])
 result=subprocess.run(['ldd',str(p)],capture_output=True,text=True,env=env)
 if 'not found' in result.stdout:raise ValueError('Unresolved native dependency '+str(p))
 for line in result.stdout.splitlines():
  m=re.search(r'=> (/\S+)',line)
  if not m:continue
  lib=Path(m[1]);leaf=lib.name
  if leaf.startswith(('libc.so','libm.so','libpthread.so','librt.so','libdl.so','ld-linux')):continue
  sha=digest(lib)
  if leaf in libs and libs[leaf]['sha256']!=sha:raise ValueError('Conflicting native library '+leaf)
  libs[leaf]=dict(source=str(lib),sha256=sha)
with tarfile.open(archive,'w') as tar:
 def add_tree(source,target):
  source=Path(source).resolve()
  if source.is_dir():vendor_dirs.extend(str(d) for d in source.rglob('*') if d.is_dir() and d.name.endswith('.libs'))
  inventory.append(dict(source=str(source),target=target))
  tar.add(source,arcname=target)
  if source.is_dir():
   for p in source.rglob('*'):native_libraries(p)
  else:native_libraries(source)
 for name in ['STAR','samtools','bedtools','featureCounts','pigz']:
  p=shutil.which(name)
  if not p:raise ValueError('Missing required module binary '+name)
  add_tree(p,'bin/'+name)
 bowtie=Path(shutil.which('bowtie2')).resolve().parent
 for p in bowtie.glob('bowtie2*'):add_tree(p,'bin/'+p.name)
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/cutadapt/4.8/lib/python3.9/site-packages','cutadapt_site')
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/cutadapt/4.8/bin/cutadapt','entrypoints/cutadapt.py')
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/MACS3/3.0.1/lib/python3.9/site-packages','macs3_site')
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/MACS3/3.0.1/bin/macs3','entrypoints/macs3.py')
 add_tree('/stornext/System/data/software/rhel/9/base/tools/openjdk/21.0.2','java')
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/picard-tools/3.1.1/picard.jar','picard/picard.jar')
 add_tree('/stornext/System/data/software/rhel/9/base/bioinf/fastqc/0.12.1','fastqc')
 add_tree('/usr/bin/perl','perl/bin/perl')
 perl_paths=subprocess.check_output(['/usr/bin/perl','-MConfig','-e','print join("\\n",map {$Config{$_}} qw(privlibexp archlibexp sitelibexp sitearchexp vendorlibexp vendorarchexp))']).decode().splitlines()
 for directory in dict.fromkeys(perl_paths):
  if Path(directory).is_dir():add_tree(directory,'perl'+directory)
 for name,data in libs.items():tar.add(Path(data['source']).resolve(),arcname='host-libs/'+name)
(root/'evidence/raw_tools_v02_inventory.json').write_text(json.dumps(dict(payloads=inventory,libraries=libs,python_runtime='Pinned separately to 3.9.18; module system Python is now 3.9.21'),indent=2)+'\n')
