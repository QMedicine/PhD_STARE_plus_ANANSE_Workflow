#!/usr/bin/env bash
#SBATCH --job-name=phd_contracts_v02
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR=${1:?Frozen source directory required}
[[ "$SOURCE_DIR" == "$TASK_ROOT"/source_* ]] || exit 2
run="$TASK_ROOT/runs/contracts_v02_$SLURM_JOB_ID";mkdir -p "$run";cd "$run"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" --env PHDSA_TEST_SOURCE="$SOURCE_DIR" \
 /vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img \
 python -m unittest discover -s "$SOURCE_DIR/tests" -v
matrix=$(awk 'NR==2 {print $2}' "$SOURCE_DIR/templates/hic_checkpoints_Old_and_Young_HSCs_v01.tsv")
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hicexplorer_v01.sif" \
 python - "$SOURCE_DIR/bin/validate_hic_resolution.py" "$matrix" <<'PY'
import importlib.util,json,sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('resolution',sys.argv[1]);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
module.validate(sys.argv[2],100000)
passed=json.loads(Path('hic_resolution.json').read_text())
try:module.validate(sys.argv[2],30000)
except ValueError:passed['wrong_resolution_rejected']=True
else:raise AssertionError('Wrong Hi-C resolution was accepted')
Path('hic_resolution_guard.json').write_text(json.dumps(passed,indent=2)+'\n');print(json.dumps(passed,indent=2))
PY
