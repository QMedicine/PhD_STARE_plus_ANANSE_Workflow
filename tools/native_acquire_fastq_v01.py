#!/usr/bin/env python3
"""Acquire one authenticated public FASTQ into the authorized scratch input cache."""
import argparse,gzip,hashlib,json,os,ssl,urllib.request
from pathlib import Path

def acquire(manifest,root,ca_file=None):
    record=json.loads(Path(manifest).read_text())
    root=Path(root).resolve()
    allowed=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow').resolve()
    if allowed not in root.parents:raise ValueError('FASTQ destination must be inside workflow scratch')
    url=record['url']
    if not any(url.startswith(prefix+'ftp.sra.ebi.ac.uk/vol1/fastq/') for prefix in ['https://','ftp://']):raise ValueError('Expected official ENA FASTQ URL')
    dest=root/Path(url).name;root.mkdir(parents=True,exist_ok=True)
    def identity(path):
        md5=hashlib.md5();sha=hashlib.sha256();size=0
        with path.open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):md5.update(block);sha.update(block);size+=len(block)
        return md5.hexdigest(),sha.hexdigest(),size
    if not dest.exists():
        partial=dest.with_suffix(dest.suffix+'.partial')
        # A failed transfer remains partial; do not treat its size as a valid input.
        context=ssl.create_default_context(cafile=ca_file)
        with urllib.request.urlopen(url,timeout=60,context=context) as source,partial.open('wb') as output:
            for block in iter(lambda:source.read(4*1024*1024),b''):output.write(block)
        md5,sha,size=identity(partial)
        if md5!=record['md5'] or size!=record['bytes']:
            partial.with_suffix(partial.suffix+'.failed_identity.json').write_text(json.dumps(dict(expected=record,observed_md5=md5,observed_sha256=sha,observed_bytes=size,status='FAILED_IDENTITY'),indent=2)+'\n')
            raise ValueError('Downloaded FASTQ differs from ENA identity')
        with gzip.open(partial,'rb') as compressed:
            for block in iter(lambda:compressed.read(4*1024*1024),b''):pass
        os.replace(partial,dest)
    md5,sha,size=identity(dest)
    if md5!=record['md5'] or size!=record['bytes']:raise ValueError('Cached FASTQ differs from ENA identity')
    receipt=dict(record,path=str(dest),observed_md5=md5,observed_sha256=sha,observed_bytes=size,status='VERIFIED',gzip_integrity='PASS')
    receipt['acquisition_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if ca_file:receipt['ca_bundle_sha256']=hashlib.sha256(Path(ca_file).read_bytes()).hexdigest()
    dest.with_suffix(dest.suffix+'.identity.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--root',required=True);p.add_argument('--ca-file');a=p.parse_args();acquire(a.manifest,a.root,a.ca_file)
