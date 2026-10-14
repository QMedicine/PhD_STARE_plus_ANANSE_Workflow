#!/usr/bin/env bash
#SBATCH --job-name=phd_state_tpm_reports
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module purge
module load nextflow/26.04.3
module swap apptainer/1.4.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR="$TASK_ROOT/source_state_tpm_reports_v44"
CONTRAST=${1:?Report contrast key required}
[[ "$CONTRAST" == 'Young_to_Old' || "$CONTRAST" == 'Old_to_Young' ]] || exit 2
RUN_DIR="$TASK_ROOT/runs/state_specific_reports_v44_$CONTRAST"
RESULTS_DIR="$TASK_ROOT/runs/state_specific_tpm_v42/results"
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/state_specific_reports_v44_$CONTRAST"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_state_reports_v44_$CONTRAST" APPTAINER_TMPDIR="$TASK_ROOT/tmp/state_specific_reports_v44_$CONTRAST"
export NXF_OPTS='-Xms512m -Xmx6g'
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
cd "$RUN_DIR"
nextflow run "$SOURCE_DIR/report_replay.nf" -c "$SOURCE_DIR/nextflow.config" -profile wehi \
 -params-file "$RUN_DIR/params.json" --trace_suffix "$SLURM_JOB_ID" --outdir "$RESULTS_DIR" \
 -work-dir "$RUN_DIR/work" -resume
