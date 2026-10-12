#!/usr/bin/env bash
#SBATCH --job-name=phd_fixed_contacts_stare
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load nextflow/26.04.3
module swap apptainer/1.4.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR="$TASK_ROOT/source_v26"
RUN_DIR="$TASK_ROOT/runs/fixed_contacts_stare_smoke_v26"
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/fixed_contacts_stare_smoke_v25"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_fixed_contacts_stare_smoke_v25" APPTAINER_TMPDIR="$TASK_ROOT/tmp/fixed_contacts_stare_smoke_v25"
export NXF_OPTS='-Xms512m -Xmx6g'
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
cd "$RUN_DIR"
nextflow run "$SOURCE_DIR/tests/native_fixed_contacts_stare_smoke.nf" -c "$SOURCE_DIR/nextflow.config" -profile wehi \
 -params-file "$TASK_ROOT/inputs/fixed_contacts_stare_fixture_v01/params.json" --outdir "$RUN_DIR/results" -work-dir "$RUN_DIR/work"
