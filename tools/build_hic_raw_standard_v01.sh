#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_raw_standard
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
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/hic_raw_standard_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/hic_raw_standard_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" containers/hic_python_wheels_v01
[[ ! -e containers/hic_raw_standard_v01.sif ]] || { echo 'Refuse image overwrite';exit 2; }
sha256sum --check evidence/raw_tools_v07.sif.sha256
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" containers/raw_tools_v07.sif python -m pip download \
 --cert "$TASK_ROOT/containers/trust/ca_bundle_v01.pem" --no-deps --only-binary=:all: \
 --dest "$TASK_ROOT/containers/hic_python_wheels_v01" biopython==1.83 pysam==0.22.1 numpy==1.26.4
sha256sum containers/hic_python_wheels_v01/*.whl > evidence/hic_python_wheels_v01.sha256
apptainer build --fakeroot containers/hic_raw_standard_v01.sif containers/hic_raw_standard_v01.def
sha256sum containers/hic_raw_standard_v01.sif > evidence/hic_raw_standard_v01.sif.sha256
apptainer inspect --json containers/hic_raw_standard_v01.sif > evidence/hic_raw_standard_v01.inspect.json
