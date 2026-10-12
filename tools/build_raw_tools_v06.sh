#!/usr/bin/env bash
#SBATCH --job-name=phd_raw_tools_v06
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
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/raw_tools_v06_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/raw_tools_v06_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
[[ ! -e containers/raw_tools_v06.sif ]] || { echo 'Refusing image overwrite';exit 2; }
sha256sum --check evidence/raw_tools_v05.sif.sha256
python3 tools/package_procps_v01.py
apptainer build --fakeroot containers/raw_tools_v06.sif containers/raw_tools_v06.def
sha256sum containers/raw_tools_v06.sif > evidence/raw_tools_v06.sif.sha256
apptainer inspect --json containers/raw_tools_v06.sif > evidence/raw_tools_v06.inspect.json
