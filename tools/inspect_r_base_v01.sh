set -euo pipefail
module load apptainer/1.3.5
SIF=/vast/scratch/users/stacpoole.q/FateFactors_development/runs/stage009-r-base-0342300-20261003-r1/r_base.sif
sha256sum "$SIF"
apptainer inspect --json "$SIF"
apptainer exec --cleanenv "$SIF" cat /etc/os-release
apptainer exec --cleanenv "$SIF" Rscript --vanilla - <<'RS'
cat(R.version.string, '\n')
p <- installed.packages(); write.table(p[,c('Package','Version')],row.names=FALSE,quote=FALSE,sep='\t')
RS
