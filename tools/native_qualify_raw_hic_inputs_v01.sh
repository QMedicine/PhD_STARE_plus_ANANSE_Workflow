#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_raw_inputs
#SBATCH --array=0-11%4
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" \
 python "$TASK_ROOT/tools/qualify_raw_hic_cohort_v01.py" --samplesheet "$TASK_ROOT/source_v22/templates/samplesheet_HiC_full_raw_cohort_v01.tsv" \
 --task-index "$SLURM_ARRAY_TASK_ID" --raw-results "$TASK_ROOT/runs/hic_full_raw_cohort_v22/results" \
 --outdir "$TASK_ROOT/evidence/raw_hic_inputs_v01" --inputs-only
