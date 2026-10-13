#!/usr/bin/env bash
#SBATCH --job-name=phd_report_publish
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$TASK_ROOT/source_v33/tools/qualify_selected_reports_v01.py" \
 --results "$TASK_ROOT/runs/selected_reports_v32/results" \
 --manifest "$PROJECT_RUN/report_manifest_selected_defaults_v01.tsv" \
 --destination "$PROJECT_RUN/validated_results/selected_reports_v32"
