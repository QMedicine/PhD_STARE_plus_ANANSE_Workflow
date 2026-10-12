#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_presplit_fixture
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" \
 python "$TASK_ROOT/tools/native_make_hic_presplit_fixture_v01.py"
