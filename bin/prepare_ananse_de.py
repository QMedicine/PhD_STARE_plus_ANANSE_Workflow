#!/usr/bin/env python3
"""Preserve DE data bytes while avoiding ANANSE 0.5.1's named-index mapping bug."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path

def prepare(source,output):
    source,output=Path(source),Path(output)
    if source.resolve()==output.resolve():raise ValueError('DE compatibility output must be separate from its input')
    with source.open(newline='') as stream:
        reader=csv.reader(stream,delimiter='\t');head=next(reader,[])
        if not head or head[0] not in ['','gene_id','gene_name','Geneid','gene','index']:raise ValueError('DE first column must identify genes')
        if len(set(head))!=len(head):raise ValueError('Duplicate DE column headings')
        seen=set()
        for row in reader:
            if len(row)!=len(head) or not row[0] or row[0] in seen:raise ValueError('Duplicate, missing or malformed DE gene identities')
            seen.add(row[0])
        if not seen:raise ValueError('Empty DE table')
    with source.open('rb') as src,output.open('wb') as dst:
        header=src.readline();fields=header.split(b'\t')
        if len(fields)<3 or b'log2FoldChange' not in fields or not any(f.strip()==b'padj' for f in fields):raise ValueError('DE columns missing')
        previous=fields[0].decode();fields[0]=b'index';dst.write(b'\t'.join(fields));shutil.copyfileobj(src,dst)
    def body_sha(p):
        h=hashlib.sha256()
        with p.open('rb') as f:
            f.readline()
            for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
        return h.hexdigest()
    input_body,output_body=body_sha(source),body_sha(output)
    if input_body!=output_body:raise ValueError('DE data bytes changed')
    Path('de_compatibility.json').write_text(json.dumps(dict(original_index_heading=previous,compatibility_index_heading='index',data_bytes_unchanged=True,data_body_sha256=input_body,reason='ANANSE gene mapping resets index then groups by the literal index column'),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args();prepare(a.input,a.output)
