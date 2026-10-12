#!/usr/bin/env python3
"""Resume an interrupted official ENA FTP transfer; accept only verified complete FASTQ."""
import argparse,ftplib,gzip,hashlib,json,os,time,urllib.parse
from pathlib import Path

def acquire(manifest,root):
    record=json.loads(Path(manifest).read_text());root=Path(root).resolve()
    if Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow') not in root.parents:raise ValueError('Unexpected scratch destination')
    parsed=urllib.parse.urlparse(record['url'])
    if parsed.scheme!='ftp' or parsed.hostname!='ftp.sra.ebi.ac.uk' or not parsed.path.startswith('/vol1/fastq/'):raise ValueError('Expected official ENA FTP path')
    root.mkdir(parents=True,exist_ok=True);dest=root/Path(parsed.path).name;partial=dest.with_suffix(dest.suffix+'.partial')
    ledger=partial.with_suffix(partial.suffix+'.attempts.json');attempts=[]
    if ledger.exists():attempts=json.loads(ledger.read_text())
    if not dest.exists():
        for attempt in range(3):
            offset=partial.stat().st_size if partial.exists() else 0
            if offset>record['bytes']:raise ValueError('Partial file is larger than the authenticated object')
            if offset==record['bytes']:break
            begin=time.monotonic();last=begin
            try:
                with ftplib.FTP(parsed.hostname,timeout=180) as ftp,partial.open('ab') as output:
                    ftp.login();ftp.voidcmd('TYPE I')
                    if ftp.size(parsed.path)!=record['bytes']:raise ValueError('ENA FTP object size changed')
                    received=offset
                    def write(block):
                        nonlocal received,last
                        output.write(block);received+=len(block)
                        if time.monotonic()-last>=60:
                            print('PROGRESS',received,'of',record['bytes'],flush=True);last=time.monotonic()
                    print('RESUME_OFFSET',offset,flush=True)
                    ftp.retrbinary('RETR '+parsed.path,write,blocksize=1024*1024,rest=offset if offset else None)
                    output.flush();os.fsync(output.fileno())
                attempts.append(dict(offset=offset,status='TRANSFER_COMPLETE',elapsed_seconds=round(time.monotonic()-begin,2)))
                ledger.write_text(json.dumps(attempts,indent=2)+'\n');break
            except (OSError,EOFError,ftplib.Error) as e:
                attempts.append(dict(offset=offset,status='INTERRUPTED',error=type(e).__name__+': '+str(e),elapsed_seconds=round(time.monotonic()-begin,2)))
                ledger.write_text(json.dumps(attempts,indent=2)+'\n');print('INTERRUPTED',attempts[-1],flush=True)
                if attempt==2:raise
                time.sleep(10*(attempt+1))
    candidate=dest if dest.exists() else partial
    md5=hashlib.md5();sha=hashlib.sha256();n=0
    with candidate.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):md5.update(block);sha.update(block);n+=len(block)
    receipt=dict(record,path=str(dest),observed_md5=md5.hexdigest(),observed_sha256=sha.hexdigest(),observed_bytes=n,attempts=attempts,acquisition_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if n!=record['bytes'] or md5.hexdigest()!=record['md5']:
        receipt['status']='FAILED_IDENTITY';partial.with_suffix(partial.suffix+'.failed_identity.json').write_text(json.dumps(receipt,indent=2)+'\n');raise ValueError('Complete file differs from ENA identity; retain for investigation')
    with gzip.open(candidate,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):pass
    if candidate==partial:os.replace(partial,dest)
    receipt.update(status='VERIFIED',gzip_integrity='PASS');dest.with_suffix(dest.suffix+'.identity.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--root',required=True);a=p.parse_args();acquire(a.manifest,a.root)
