#!/usr/bin/env bash
#SBATCH --job-name=phd_source_snapshot
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
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" python - <<'PY'
from pathlib import Path
import hashlib,json,shutil,tarfile
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');source=root/'source_v12'
scratch=root/'evidence/source_snapshot_v12';scratch.mkdir(parents=True,exist_ok=True)
files=sorted(p for p in source.rglob('*') if p.is_file())
if not files or sum(p.stat().st_size for p in files)>20*1024*1024:raise ValueError('Unexpected source snapshot size')
archive=scratch/'source_v12.tar.gz'
with tarfile.open(archive,'w:gz') as out:
 for p in files:out.add(p,arcname=str(Path('source_v12')/p.relative_to(source)),recursive=False)
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
inventory=scratch/'source_v12.sha256';inventory.write_text('\n'.join(digest(p)+'  '+str(p.relative_to(source)) for p in files)+'\n')
target=Path('/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/validated_results/atac_saved_union_v12/provenance')
for p in [archive,inventory]:
 dest=target/p.name
 if dest.exists() and digest(dest)!=digest(p):raise ValueError('Refuse to replace differing source snapshot')
 if not dest.exists():shutil.copy2(p,dest)
 if digest(dest)!=digest(p):raise ValueError('Source snapshot copy differs')
record=dict(source=str(source),source_files=len(files),source_bytes=sum(p.stat().st_size for p in files),archive_sha256=digest(archive),inventory_sha256=digest(inventory),durable_directory=str(target),status='SOURCE_SNAPSHOT_COPY_VERIFIED')
(scratch/'source_snapshot_receipt.json').write_text(json.dumps(record,indent=2)+'\n');shutil.copy2(scratch/'source_snapshot_receipt.json',target/'source_snapshot_receipt.json');print(json.dumps(record,indent=2))
PY
