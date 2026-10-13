#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_single_publish
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$TASK_ROOT/source_v38/tools/qualify_hic_single_library_results_v01.py" \
 --results "$TASK_ROOT/runs/hic_single_library_to_stare_v38/results" \
 --samplesheet "$PROJECT_RUN/samplesheet_HiC_single_library_practice_v01.tsv" \
 --bam-qualifications "$TASK_ROOT/evidence/hic_single_library_bams_v38" \
 --destination "$PROJECT_RUN/validated_results/hic_single_library_to_stare_v38"
