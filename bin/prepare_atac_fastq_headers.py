#!/usr/bin/env python3
"""Explicit SRA spot/length header reconstruction; sequences and qualities stay intact."""
import argparse,gzip,hashlib,json,re
from pathlib import Path

def prepare(paths,prefix):
    receipts=[]
    for mate,path in enumerate(paths,1):
        path=Path(path);output=Path(f'{prefix}.R{mate}.fastq');digest=hashlib.sha256();records=0
        open_input=gzip.open if path.name.endswith('.gz') else open
        with open_input(path,'rb') as source,output.open('wb') as target:
            while True:
                header=source.readline()
                if not header:break
                sequence=source.readline();plus=source.readline();quality=source.readline()
                name=header.rstrip(b'\r\n').split(None,1)[0]
                match=re.fullmatch(rb'@((?:SRR|ERR|DRR)[0-9]+\.([0-9]+))(?:/[12])?',name)
                if not match or not plus.startswith(b'+') or not quality or len(sequence.rstrip(b'\r\n'))!=len(quality.rstrip(b'\r\n')):raise ValueError('Invalid SRA-name FASTQ record')
                new_header=b'@'+match[1]+b' '+match[2]+b' length='+str(len(sequence.rstrip(b'\r\n'))).encode()+b'\n'
                block=new_header+sequence+b'+\n'+quality
                target.write(block);digest.update(block);records+=1
        receipts.append(dict(input=str(path),output=str(output),records=records,output_sha256=digest.hexdigest(),sequence_and_quality_lines_unchanged=True))
    if len(receipts)==2 and receipts[0]['records']!=receipts[1]['records']:raise ValueError('Mate record counts differ')
    Path('header_transformation.json').write_text(json.dumps(dict(method='sra_spot_length',format='@ACCESSION.SPOT SPOT length=READ_LENGTH',reads=receipts,historical_equivalence='REQUIRES_COMPLETE_ALIGNMENT_REGRESSION'),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prefix',required=True);p.add_argument('--reads',nargs='+',required=True);a=p.parse_args();prepare(a.reads,a.prefix)
