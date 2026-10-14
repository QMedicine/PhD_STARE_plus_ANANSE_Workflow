#!/usr/bin/env bash
#SBATCH --job-name=phd_state_tpm_qualify
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module purge
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR="$TASK_ROOT/source_state_tpm_reports_v44"
RUN_DIR="$TASK_ROOT/runs/state_specific_tpm_v42"
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$SOURCE_DIR/tools/qualify_state_tpm_run_v01.py" --results "$RUN_DIR/results" \
 --params "$RUN_DIR/params.json" --input-receipt "$RUN_DIR/input_receipt.json" \
 --destination "$PROJECT_RUN/validated_results/state_specific_tpm_v42"
