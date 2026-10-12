#!/usr/bin/env bash
#SBATCH --job-name=phd_contracts_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/runs/contracts_v01_$SLURM_JOB_ID";mkdir -p "$run";cd "$run"
export APPTAINERENV_PHDSA_TEST_SOURCE="$TASK_ROOT/contract_validation_source_v01"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" /vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img \
 python -m unittest discover -s "$TASK_ROOT/contract_validation_source_v01/tests" -v
