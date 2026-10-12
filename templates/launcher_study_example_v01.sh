#!/usr/bin/env bash
#SBATCH --job-name=phd_stare_plus_ananse
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
trap 'printf "ERROR line=%s command=%q\n" "$LINENO" "$BASH_COMMAND" >&2' ERR
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
SOURCE_DIR=${PHDSA_SOURCE_DIR:?Set PHDSA_SOURCE_DIR to a frozen, validated scratch source checkout}
PARAMS_FILE=${PHDSA_PARAMS_FILE:?Set PHDSA_PARAMS_FILE to the reviewed run JSON}
[[ "$SOURCE_DIR" == "$TASK_ROOT"/source_* ]] || { echo 'Source must be a versioned scratch checkout'; exit 2; }
[[ "$PARAMS_FILE" == "$PROJECT_RUN"/* || "$PARAMS_FILE" == "$TASK_ROOT"/* ]] || { echo 'Unexpected parameter file location'; exit 2; }
[[ -r "$SOURCE_DIR/main.nf" && -r "$PARAMS_FILE" ]] || { echo 'Source or parameters missing'; exit 2; }
module load nextflow/26.04.3 apptainer/1.3.5
RUN_ID=${PHDSA_RUN_ID:-Old_and_Young_HSCs_v01_${SLURM_JOB_ID:-manual}}
[[ "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid run ID'; exit 2; }
RUN_DIR="$TASK_ROOT/runs/$RUN_ID"
mkdir -p "$RUN_DIR/provenance" "$TASK_ROOT/cache" "$TASK_ROOT/tmp/$RUN_ID"
export NXF_HOME="$TASK_ROOT/cache/nextflow"
export NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_$RUN_ID"
export APPTAINER_TMPDIR="$TASK_ROOT/tmp/$RUN_ID"
export NXF_OPTS='-Xms512m -Xmx6g'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RUN_DIR"
nextflow -version > provenance/nextflow_version.txt
apptainer --version > provenance/apptainer_version.txt
sha256sum "$PARAMS_FILE" > provenance/params.sha256
cp "$PARAMS_FILE" provenance/params.json
nextflow config "$SOURCE_DIR" -profile wehi -o flat > provenance/profile_defaults.config
cmd=(nextflow run "$SOURCE_DIR/main.nf" -profile wehi -params-file "$PARAMS_FILE" -work-dir "$RUN_DIR/work" --outdir "$RUN_DIR/results" --trace_suffix "${SLURM_JOB_ID:-manual}")
if [[ ${PHDSA_RESUME:-FALSE} == TRUE ]]; then cmd+=(-resume); fi
printf 'COMMAND '; printf '%q ' "${cmd[@]}"; printf '\n'
"${cmd[@]}"
# The final copy is performed by the separately validated publication helper.
# A successful driver alone is not an exact-reproduction receipt.
