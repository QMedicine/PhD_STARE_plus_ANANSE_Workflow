#!/usr/bin/env bash
#SBATCH --job-name=phd_atac_paired_fixture
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/raw_tools_v06.sif" \
 python "$TASK_ROOT/tools/native_make_paired_atac_fixture_v01.py" --source "$TASK_ROOT/source_v14" --root "$TASK_ROOT/inputs/paired_atac_fixture_v01"
