#!/usr/bin/env bash
#SBATCH --job-name=phd_r_container_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load R/4.5.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/build_r_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/build_r_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" evidence
[[ ! -e containers/r_plot_v01.sif ]] || { echo 'Refusing image overwrite'; exit 2; }
python3 tools/package_r_runtime_v01.py
sha256sum containers/r_runtime_v01.tar > evidence/r_runtime_v01.tar.sha256
apptainer build --fakeroot containers/r_plot_v01.sif containers/r_plot_v01.def
sha256sum containers/r_plot_v01.sif > evidence/r_plot_v01.sif.sha256
apptainer inspect --json containers/r_plot_v01.sif > evidence/r_plot_v01.inspect.json
printf 'BUILD_COMPLETE\n'
