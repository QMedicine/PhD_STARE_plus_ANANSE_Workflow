#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_single_library
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
export PHDSA_PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
export PHDSA_SOURCE_DIR=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/source_v38
export PHDSA_PARAMS_FILE="$PHDSA_PROJECT_RUN/params_HiC_single_library_to_STARE_v01.json"
export PHDSA_RUN_ID=hic_single_library_to_stare_v38
exec bash "$PHDSA_SOURCE_DIR/templates/launcher_template_v03.sh"
