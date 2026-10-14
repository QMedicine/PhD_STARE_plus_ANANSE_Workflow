#!/usr/bin/env bash
#SBATCH --job-name=phd_state_tpm_receipt
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module purge
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "${1:-$TASK_ROOT/transfer/finalize_state_tpm_receipt_v01.py}" \
 --results "$TASK_ROOT/runs/state_specific_tpm_v42/results" \
 --output "$TASK_ROOT/runs/state_specific_tpm_v42/publication_receipt_check.json"
