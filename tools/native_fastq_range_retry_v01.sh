#!/usr/bin/env bash
#SBATCH --job-name=phd_fastq_range_retry
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/raw_tools_v06.sif" python - <<'PY'
import gzip,hashlib,json,os,shutil,ssl,urllib.request
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
manifest=json.loads((root/'source_v14/templates/ena_ATAC_pilot_identity_v01.json').read_text())
diagnostic=json.loads((root/'evidence/fastq_byte_diagnostic_v01/byte_diagnostic.json').read_text())
folder=root/'inputs/public_fastq_range_v03';folder.mkdir(parents=True,exist_ok=True)
candidate=folder/'SRR13186609.fastq.gz.partial'
if candidate.exists():raise ValueError('Range candidate already exists; preserve it')
shutil.copyfile(diagnostic['paths'][1],candidate)
context=ssl.create_default_context(cafile=str(root/'containers/trust/ca_bundle_v01.pem'))
receipts=[]
for chunk in diagnostic['different_chunks']:
 start,end=chunk['start'],chunk['end']-1
 request=urllib.request.Request(manifest['url'],headers={'Range':f'bytes={start}-{end}'})
 with urllib.request.urlopen(request,timeout=180,context=context) as response:
  content_range=response.headers.get('Content-Range');payload=response.read(end-start+2)
  if response.status!=206 or content_range!=f'bytes {start}-{end}/{manifest["bytes"]}' or len(payload)!=end-start+1:raise ValueError('Server did not return the exact requested byte range')
 with candidate.open('r+b') as output:output.seek(start);output.write(payload);output.flush();os.fsync(output.fileno())
 receipts.append(dict(start=start,end=end,bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),content_range=content_range))
 print('RANGE_RETRIEVED',start,end,flush=True)
md5=hashlib.md5();sha=hashlib.sha256()
with candidate.open('rb') as f:
 for block in iter(lambda:f.read(4*1024*1024),b''):md5.update(block);sha.update(block)
record=dict(manifest,observed_md5=md5.hexdigest(),observed_sha256=sha.hexdigest(),observed_bytes=candidate.stat().st_size,ranges=receipts,method='Official HTTP Range retrieval into a separate candidate; global ENA checksum remains mandatory')
if record['observed_md5']!=manifest['md5'] or record['observed_bytes']!=manifest['bytes']:
 record['status']='FAILED_IDENTITY';(folder/'range_retry.json').write_text(json.dumps(record,indent=2)+'\n');raise ValueError('Range retry still differs from ENA MD5; preserve candidate')
with gzip.open(candidate,'rb') as f:
 for block in iter(lambda:f.read(4*1024*1024),b''):pass
dest=folder/'SRR13186609.fastq.gz';os.replace(candidate,dest);record.update(status='VERIFIED',gzip_integrity='PASS',path=str(dest));dest.with_suffix(dest.suffix+'.identity.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
PY
