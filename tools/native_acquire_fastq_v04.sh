#!/usr/bin/env bash
#SBATCH --job-name=phd_fastq_resume_ftp
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/raw_tools_v06.sif" \
 python "$TASK_ROOT/tools/native_acquire_fastq_v04.py" --manifest "$1" --root "$TASK_ROOT/inputs/public_fastq_ftp_v02"
