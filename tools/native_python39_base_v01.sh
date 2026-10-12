#!/usr/bin/env bash
#SBATCH --job-name=phd_python39_base
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/runs/python39_base_v01_$SLURM_JOB_ID";mkdir -p "$run";cd "$run"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/python39_base_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/python39_base_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
image="$TASK_ROOT/containers/python_3.9.18_slim_bookworm_v01.sif"
[[ ! -e "$image" ]] || { echo 'Refusing image overwrite';exit 2; }
apptainer pull --arch amd64 "$image" docker://python:3.9.18-slim-bookworm
apptainer exec --cleanenv --no-home "$image" python -c 'import sys,platform;assert sys.version_info[:3]==(3,9,18);assert platform.machine()=="x86_64";print(sys.version)'
sha256sum "$image" > "$TASK_ROOT/evidence/python39_base_v01.sif.sha256"
apptainer inspect --json "$image" > "$TASK_ROOT/evidence/python39_base_v01.inspect.json"
printf 'PYTHON39_BASE_COMPLETE\n'
