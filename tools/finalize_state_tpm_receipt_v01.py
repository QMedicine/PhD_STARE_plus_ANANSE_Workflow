#!/usr/bin/env python3
"""Verify published bytes and exclude a finalized receipt's stale self-hash."""
import argparse
import hashlib
import json
from pathlib import Path

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''):digest.update(block)
    return digest.hexdigest()

def finalize(results,output):
    root=Path(results);receipt=root/'qualification.json'
    record=json.loads(receipt.read_text())
    if record['publication_status']!='CHECKSUM_VERIFIED_COPY':
        raise ValueError('Published release is not qualified')
    destination=Path(record['published_destination'])
    if sha(receipt)!=sha(destination/'qualification.json'):
        raise ValueError('Source and destination receipts differ')
    checks=record['published_file_sha256']
    checks.pop('qualification.json',None)
    for relative,digest in checks.items():
        if sha(root/relative)!=digest or sha(destination/relative)!=digest:
            raise ValueError('Published file differs: '+relative)
    record['receipt_hash_policy']='Finalized qualification receipt excluded from its own payload hash map; external receipt checksum records its final bytes'
    text=json.dumps(record,indent=2)+'\n'
    receipt.write_text(text);(destination/'qualification.json').write_text(text)
    assert sha(receipt)==sha(destination/'qualification.json')
    verification=dict(status='VERIFIED_PUBLISHED_STATE_SPECIFIC_TPM_RELEASE',
        payload_files=len(checks),receipt_sha256=sha(receipt),published_destination=str(destination),
        source_receipt=str(receipt),self_reference_excluded=True)
    Path(output).write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--results',required=True);parser.add_argument('--output',required=True)
    finalize(**vars(parser.parse_args()))
