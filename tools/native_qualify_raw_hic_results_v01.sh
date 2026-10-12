#!/usr/bin/env bash
#SBATCH --job-name=phd_raw_hic_publish
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" python "$TASK_ROOT/tools/qualify_raw_hic_results_v01.py"
