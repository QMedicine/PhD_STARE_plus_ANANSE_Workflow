#!/usr/bin/env bash
#SBATCH --job-name=phd_state_specific_tpm
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
RUN_DIR="$TASK_ROOT/runs/state_specific_tpm_v42"
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/state_specific_tpm_v42"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_state_specific_tpm_v42" APPTAINER_TMPDIR="$TASK_ROOT/tmp/state_specific_tpm_v42"
export NXF_OPTS='-Xms512m -Xmx6g'
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
printf 'job=%s host=%s source=%s run=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$SOURCE_DIR" "$RUN_DIR"
date -Is
cd "$RUN_DIR"
nextflow run "$SOURCE_DIR/main.nf" -c "$SOURCE_DIR/nextflow.config" -profile wehi \
  -params-file "$RUN_DIR/params.json" --trace_suffix "$SLURM_JOB_ID" \
  --outdir "$RUN_DIR/results" -work-dir "$RUN_DIR/work" -resume
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
  python "$SOURCE_DIR/tools/qualify_state_tpm_run_v01.py" --results "$RUN_DIR/results" \
  --params "$RUN_DIR/params.json" --input-receipt "$RUN_DIR/input_receipt.json" \
  --destination "$PROJECT_RUN/validated_results/state_specific_tpm_v42"
