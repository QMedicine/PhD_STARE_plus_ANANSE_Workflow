set -euo pipefail
module load miniconda3
printf '\n=== HOST ===\n'; hostname; uname -m; date -u
printf '\n=== MODULES ===\n'; module avail nextflow apptainer R STAR subread bowtie2 samtools bedtools macs3 2>&1
for env_name in ananse cooler-env hicexplorer python3; do
 printf '\n=== CONDA %s ===\n' "$env_name"
 conda list -n "$env_name" --json
 printf '\n=== EXPLICIT %s ===\n' "$env_name"
 conda list -n "$env_name" --explicit
 done
printf '\n=== R PACKAGES ===\n'
module load R/4.5.1
R_LIBS_USER=/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1 Rscript --vanilla -e 'cat(R.version.string,"\n"); x<-installed.packages(); write.table(x[,c("Package","Version","LibPath","Built")],row.names=FALSE,sep="\t",quote=FALSE)'
printf '\n=== SLURM ===\n'; sinfo -o '%P %l %c %m'; sacctmgr -n show assoc where user=stacpoole.q format=Account,Partition,QOS
printf '\n=== TF ACTIVITY REVISION ===\n'; git -C /vast/projects/Immune_Aging/.nextflow_assets/nf-core/tfactivity rev-parse HEAD; git -C /vast/projects/Immune_Aging/.nextflow_assets/nf-core/tfactivity status --short
