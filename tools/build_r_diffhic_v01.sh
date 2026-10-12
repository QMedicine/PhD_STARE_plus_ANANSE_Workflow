#!/usr/bin/env bash
#SBATCH --job-name=phd_r_diffhic_build
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load R/4.5.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/build_diffhic_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/build_diffhic_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
[[ ! -e containers/r_diffhic_v01.sif && ! -e containers/r_diffhic_runtime_v01.tar ]] || { echo 'Refuse overwrite';exit 2; }
Rscript --vanilla tools/remote_diffhic_library_closure_v01.R > containers/r_diffhic_library_closure_v01.tsv
python3 tools/package_diffhic_runtime_v01.py
sha256sum --check evidence/r_plot_v01.sif.sha256
apptainer build --fakeroot containers/r_diffhic_v01.sif containers/r_diffhic_v01.def
sha256sum containers/r_diffhic_v01.sif > evidence/r_diffhic_v01.sif.sha256
