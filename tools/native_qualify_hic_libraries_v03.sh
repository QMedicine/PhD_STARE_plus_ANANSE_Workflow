#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_library_validate
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR=${PHDSA_SOURCE_DIR:?Set a frozen source directory}
SAMPLESHEET=${PHDSA_SAMPLESHEET:?Set the original paired FASTQ samplesheet}
RAW_RESULTS=${PHDSA_RAW_RESULTS:?Set the raw alignment result directory}
QUALIFICATION_DIR=${PHDSA_QUALIFICATION_DIR:?Set the scratch qualification directory}
PHASE=${PHDSA_PHASE:?Set inputs or bam}
[[ "$SOURCE_DIR" == "$TASK_ROOT"/source_* && "$RAW_RESULTS" == "$TASK_ROOT"/runs/* && "$QUALIFICATION_DIR" == "$TASK_ROOT"/evidence/* ]] || exit 2
args=(--samplesheet "$SAMPLESHEET" --task-index "${SLURM_ARRAY_TASK_ID:?Array index required}" --raw-results "$RAW_RESULTS" --outdir "$QUALIFICATION_DIR")
case "$PHASE" in
 inputs) args+=(--inputs-only) ;;
 bam) args+=(--input-qualification-root "${PHDSA_INPUT_QUALIFICATION_DIR:?Set input qualification directory}") ;;
 *) echo 'Unsupported validation phase' >&2; exit 2 ;;
esac
printf 'SOURCE=%s SAMPLESHEET=%s PHASE=%s TASK=%s OUTPUT=%s\n' "$SOURCE_DIR" "$SAMPLESHEET" "$PHASE" "$SLURM_ARRAY_TASK_ID" "$QUALIFICATION_DIR"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" \
 python "$SOURCE_DIR/tools/qualify_raw_hic_cohort_v02.py" "${args[@]}"
