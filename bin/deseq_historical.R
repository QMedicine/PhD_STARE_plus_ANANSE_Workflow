#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DESeq2);library(BiocParallel)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==8)
counts<-read.delim(args[1],row.names=1,check.names=FALSE)
samples<-read.delim(args[2],check.names=FALSE)
samples<-samples[samples$assay=='RNA',]
stopifnot(!anyDuplicated(samples$sample_name),setequal(colnames(counts),samples$sample_name))
# Preserve manifest order, as in the archived differential-expression template.
counts<-counts[,samples$sample_name,drop=FALSE]
state_order<-if(nzchar(args[8]))strsplit(args[8],',',fixed=TRUE)[[1]] else unique(samples$cell_state)
stopifnot(!anyDuplicated(state_order),setequal(state_order,unique(samples$cell_state)))
samples$contrast_condition<-factor(samples$cell_state,levels=state_order)
rownames(samples)<-samples$sample_name
source<-args[3];target<-args[4]
stopifnot(source %in% state_order,target %in% state_order,source!=target)
dds<-DESeqDataSetFromMatrix(countData=round(as.matrix(counts)),colData=samples,design=~0+contrast_condition)
dds<-DESeq(dds,test='Wald',fitType='parametric',minReplicatesForReplace=7,useT=FALSE,sfType='ratio',parallel=TRUE,BPPARAM=MulticoreParam(as.integer(args[7])))
res<-results(dds,lfcThreshold=0,altHypothesis='greaterAbs',independentFiltering=TRUE,alpha=.1,pAdjustMethod='BH',minmu=.5,contrast=c('contrast_condition',as.character(target),as.character(source)))
# The archived template replaces its first results object with ashr lfcShrink.
res<-lfcShrink(dds,type='ashr',contrast=c('contrast_condition',as.character(target),as.character(source)))
out<-data.frame(gene_id=rownames(res),as.data.frame(res),check.names=FALSE)
write.table(out,args[5],sep='\t',quote=FALSE,row.names=FALSE)
reverse<-out;reverse$log2FoldChange<- -reverse$log2FoldChange
write.table(reverse,args[6],sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(R.version.string,paste('DESeq2',packageVersion('DESeq2')),paste('ashr',packageVersion('ashr'))),'versions.txt')
