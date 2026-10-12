#!/usr/bin/env bash
#SBATCH --job-name=phd_cached_container_audit_v01
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
cd "$TASK_ROOT"
mkdir -p evidence/cached_containers
for image in /vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod*data.img; do
 name=$(basename "$image")
 apptainer inspect --json "$image" > "evidence/cached_containers/${name}.inspect.json"
 if apptainer exec --cleanenv "$image" python3 -c 'import pandas,sys; print(sys.version); print(pandas.__version__)' > "evidence/cached_containers/${name}.python.txt" 2>&1; then
  sha256sum "$image" > "evidence/cached_containers/${name}.sha256"
 fi
done
image=/vast/projects/Immune_Aging/apptainer_singularity_cache/depot.galaxyproject.org-singularity-bioconductor-deseq2-1.34.0--r41hc247a5b_3.img
sha256sum "$image" > evidence/cached_containers/deseq2.sha256
apptainer exec --cleanenv "$image" Rscript --vanilla -e 'cat(R.version.string,"\n");for(p in c("DESeq2","ashr","BiocParallel")){library(p,character.only=TRUE);cat(p,as.character(packageVersion(p)),"\n")}' > evidence/cached_containers/deseq2.versions.txt
printf 'AUDIT_COMPLETE\n'
