#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(SummarizedExperiment);library(csaw)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
C <- fread(args[1]);base <- fread(args[2],col.names=c('#chrom','start','end','peak_id'))
stopifnot(identical(C$peak_id,base$peak_id),!anyDuplicated(C$peak_id),!anyDuplicated(names(C)),
          ncol(C)>1,!any(names(C)[-1] %in% names(base)))
mat <- as.matrix(C[,-1,with=FALSE]);storage.mode(mat)<-'integer'
stopifnot(!anyNA(mat),all(mat>=0),all(colSums(mat)>0))
se <- SummarizedExperiment(assays=list(counts=mat));colData(se)$totals<-colSums(mat)
se <- normOffsets(se)
cpm_mat <- calculateCPM(se,use.offsets=TRUE,use.norm.factors=FALSE,log=FALSE,assay.id='counts')
res <- cbind(base,as.data.table(cpm_mat))
fwrite(res,args[3],sep='\t',quote=FALSE)
off<-assay(se,'offset')
fwrite(data.table(sample=colnames(mat),lib.size=colData(se)$totals,offset_median=apply(off,2,median),offset_IQR=apply(off,2,function(x)diff(quantile(x,c(.25,.75))))),'loessOffsets_qc.tsv',sep='\t')
writeLines(c(R.version.string,paste('csaw',packageVersion('csaw')),paste('SummarizedExperiment',packageVersion('SummarizedExperiment')),paste('data.table',packageVersion('data.table'))),'versions.txt')
