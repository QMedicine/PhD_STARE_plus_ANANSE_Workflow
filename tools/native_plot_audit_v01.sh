#!/usr/bin/env bash
#SBATCH --job-name=phd_plot_audit
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
AUDIT_ROOT="$TASK_ROOT/runs/plot_audit_v01"
python3 "$TASK_ROOT/tools/make_plot_audit_v01.py" --source "$TASK_ROOT/source_v22" --out "$AUDIT_ROOT"
python3 "$TASK_ROOT/tools/run_plot_audit_v01.py" --root "$AUDIT_ROOT" --image "$TASK_ROOT/containers/r_plot_v01.sif"
