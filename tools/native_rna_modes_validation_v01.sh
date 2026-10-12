#!/usr/bin/env bash
#SBATCH --job-name=phd_rna_modes_validate
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
run="$TASK_ROOT/evidence/rna_modes_validation_v01";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" python - <<'PY'
import hashlib,json
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/runs')
files=['recount/selected_counts.csv','tpm/counts.tpm.tsv','tpm/counts.counts_filtered.tsv','de/Old_HSC_to_Young_HSC/forward.deseq2.tsv']
checks=[]
for relative in files:
 a=root/'rna_recount_v08/results/rna'/relative;b=root/'rna_legacy_guard_v16/results/rna'/relative
 digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 checks.append(dict(path=relative,expected_sha256=digest(a),observed_sha256=digest(b),byte_identical=digest(a)==digest(b)))
if not all(r['byte_identical'] for r in checks):raise ValueError('Historical RNA route changed')
annotation=root/'rna_annotation_mode_v16/results/rna/annotation_symbol_sum'
provenance=json.loads((annotation/'annotation_counts_provenance.json').read_text())
if provenance['count_conservation']!='EXACT' or provenance['historical_equivalence']!='NOT_CLAIMED':raise ValueError('Annotation mode lacks explicit count/provenance closure')
record=dict(historical_mode=checks,annotation_mode=provenance,status='PASS_DISTINCT_METHODS_RETAINED')
Path('rna_modes_validation.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
PY
