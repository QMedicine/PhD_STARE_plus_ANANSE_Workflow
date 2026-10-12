#!/usr/bin/env bash
#SBATCH --job-name=phd_fastq_identity
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/failed_fastq_identity_v01";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/raw_tools_v06.sif" python - <<'PY'
import gzip,hashlib,json,ssl,urllib.request
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
p=root/'inputs/public_fastq/SRR13186609.fastq.gz.partial'
h=hashlib.md5();sha=hashlib.sha256();n=0
with p.open('rb') as f:
 for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b);sha.update(b);n+=len(b)
context=ssl.create_default_context(cafile=str(root/'containers/trust/ca_bundle_v01.pem'))
manifest=json.loads((root/'source_v14/templates/ena_ATAC_pilot_identity_v01.json').read_text())
fresh=urllib.request.urlopen(manifest['source'],context=context,timeout=30).read().decode()
decoded=0;integrity=True;reason=''
try:
 with gzip.open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):decoded+=len(b)
except Exception as e:integrity=False;reason=str(e)
with urllib.request.urlopen(urllib.request.Request(manifest['url'],method='HEAD'),context=context,timeout=30) as response:
 headers={k:v for k,v in response.headers.items() if k.lower() in ['content-length','content-type','content-encoding','last-modified','etag','content-md5']}
record=dict(path=str(p),bytes=n,md5=h.hexdigest(),sha256=sha.hexdigest(),expected=manifest,fresh_ena_metadata=fresh,gzip_integrity=integrity,gzip_error=reason,decoded_bytes=decoded,http_headers=headers)
Path('failed_identity.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
PY
