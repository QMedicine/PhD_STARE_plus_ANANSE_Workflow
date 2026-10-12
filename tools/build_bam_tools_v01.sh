#!/usr/bin/env bash
#SBATCH --job-name=phd_bam_tools_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load samtools/1.22.1 bedtools/2.31.1 subread/2.0.6 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/build_bam_$SLURM_JOB_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/build_bam_$SLURM_JOB_ID"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" evidence
[[ ! -e containers/bam_tools_v01.sif ]] || { echo 'Refusing overwrite'; exit 2; }
samtools --version;bedtools --version;featureCounts -v
python3 tools/package_bam_tools_v01.py
sha256sum containers/bam_tools_v01.tar > evidence/bam_tools_v01.tar.sha256
apptainer build --fakeroot containers/bam_tools_v01.sif containers/bam_tools_v01.def
sha256sum containers/bam_tools_v01.sif > evidence/bam_tools_v01.sif.sha256
apptainer inspect --json containers/bam_tools_v01.sif > evidence/bam_tools_v01.inspect.json
printf 'BUILD_COMPLETE\n'
