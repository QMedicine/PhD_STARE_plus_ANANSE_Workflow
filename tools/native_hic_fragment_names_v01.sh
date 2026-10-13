#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_fragment_names
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
PAIR_PARAM=$(find "$TASK_ROOT/runs/hic_single_library_pilot_v36/work" -name pair_param.rds -type f -print -quit)
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/r_diffhic_v01.sif" Rscript - "$PAIR_PARAM" <<'RSCRIPT'
suppressPackageStartupMessages(library(diffHic))
x<-readRDS(commandArgs(TRUE)[1]);print(head(seqlevels(x$fragments),4));print(head(seqlengths(x$fragments),4));print(length(seqlevels(x$fragments)))
RSCRIPT
