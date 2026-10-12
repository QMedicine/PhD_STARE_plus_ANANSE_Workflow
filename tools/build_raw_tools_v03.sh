#!/usr/bin/env bash
#SBATCH --job-name=phd_raw_tools_build_v03
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load STAR/2.7.11b bowtie2/2.5.3 samtools/1.22.1 bedtools/2.31.1 subread/2.0.6 pigz apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/raw_tools_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/raw_tools_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
[[ ! -e containers/raw_tools_v03.sif ]] || { echo 'Refusing image overwrite';exit 2; }
sha256sum --check evidence/python39_base_v01.sif.sha256
python3 tools/package_raw_tools_v03.py
sha256sum containers/raw_tools_v03.tar > evidence/raw_tools_v03.tar.sha256
apptainer build --fakeroot containers/raw_tools_v03.sif containers/raw_tools_v03.def
sha256sum containers/raw_tools_v03.sif > evidence/raw_tools_v03.sif.sha256
apptainer inspect --json containers/raw_tools_v03.sif > evidence/raw_tools_v03.inspect.json
printf 'RAW_TOOLS_BUILD_COMPLETE\n'
