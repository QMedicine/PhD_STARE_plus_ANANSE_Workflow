#!/usr/bin/env bash
#SBATCH --job-name=phd_ananse_container_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/build_ananse_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/build_ananse_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" containers evidence
image="$TASK_ROOT/containers/ananse_v01.sif"
[[ ! -e "$image" ]] || { echo 'Refusing to overwrite existing image' >&2; exit 2; }
uname -m; apptainer --version
sha256sum containers/ananse_v01.def
apptainer build --fakeroot "$image" containers/ananse_v01.def
sha256sum "$image" > evidence/ananse_v01.sif.sha256
apptainer inspect --json "$image" > evidence/ananse_v01.inspect.json
apptainer exec --cleanenv --no-home "$image" ananse --version
printf 'BUILD_COMPLETE\n'
