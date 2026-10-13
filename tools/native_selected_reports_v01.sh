#!/usr/bin/env bash
#SBATCH --job-name=phd_selected_reports
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load nextflow/26.04.3
module swap apptainer/1.4.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR="$TASK_ROOT/source_v32"
RUN_DIR="$TASK_ROOT/runs/selected_reports_v32"
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/selected_reports_v32"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_selected_reports_v32" APPTAINER_TMPDIR="$TASK_ROOT/tmp/selected_reports_v32"
export NXF_OPTS='-Xms512m -Xmx6g'
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
cd "$RUN_DIR"
nextflow run "$SOURCE_DIR/report_replay.nf" -c "$SOURCE_DIR/nextflow.config" -profile wehi \
 -params-file "$SOURCE_DIR/templates/params_report_selected_defaults_v01.json" --trace_suffix "$SLURM_JOB_ID" --outdir "$RUN_DIR/results" -work-dir "$RUN_DIR/work" "$@"
