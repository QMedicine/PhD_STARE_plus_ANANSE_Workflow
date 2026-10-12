#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(diffHic);library(csaw);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<1)stop('Missing operation')
operation<-args[1]
if(operation=='fragments'){
    stopifnot(length(args)==5)
    fragments<-cutGenome(args[2],args[3],as.integer(args[4]))
    saveRDS(pairParam(fragments),args[5])
} else if(operation=='prepare'){
    stopifnot(length(args)==10)
    param<-readRDS(args[3]);prefix<-args[4]
    diagnostics<-preparePairs(args[2],param,file=paste0(prefix,'.pairs.h5'),dedup=TRUE,minq=as.integer(args[5]),chim.dist=as.integer(args[6]))
    pruning<-prunePairs(paste0(prefix,'.pairs.h5'),param,file.out=paste0(prefix,'.trimmed.h5'),max.frag=as.integer(args[7]),min.inward=as.integer(args[8]),min.outward=as.integer(args[9]))
    write_json(list(sample=args[10],diagnostics=diagnostics,pruning=pruning,parameters=list(dedup=TRUE,minq=as.integer(args[5]),chim.dist=as.integer(args[6]),max.frag=as.integer(args[7]),min.inward=as.integer(args[8]),min.outward=as.integer(args[9]))),paste0(prefix,'.qc.json'),pretty=TRUE,auto_unbox=TRUE)
} else if(operation=='count'){
    stopifnot(length(args)>=7)
    param<-readRDS(args[2]);sheet<-read.delim(args[3],stringsAsFactors=FALSE,check.names=FALSE)
    widths<-as.integer(strsplit(args[4],',',fixed=TRUE)[[1]])
    stopifnot(all(is.finite(widths)),all(widths>0),!anyDuplicated(widths))
    method<-args[5];filter<-as.integer(args[6]);files<-args[-seq_len(6)]
    labels<-sub('\\.trimmed\\.h5$','',basename(files));order<-match(sheet$sample_name[sheet$assay=='HiC'],labels)
    if(anyNA(order)||length(order)!=length(files)||anyDuplicated(labels))stop('Hi-C sample/index identity mismatch')
    files<-files[order];labels<-labels[order]
    if(length(files)<2)stop('Joint loess normalization requires at least two libraries')
    qc<-list()
    for(width in widths){
        x<-squareCounts(files,param,width=width,filter=filter)
        colnames(x)<-labels
        if(nrow(x)==0)stop('No retained interactions at resolution ',width)
        totals<-as.numeric(x$totals)
        if(any(!is.finite(totals))||any(totals<=0))stop('A library has no retained read pairs')
        if(method!='loess_offsets')stop('Unsupported diffHic normalization')
        x<-csaw::normOffsets(x)
        offsets<-assay(x,'offset')
        if(any(!is.finite(offsets)))stop('Non-finite loess offsets')
        saveRDS(x,paste0('interactions_',width,'.rds'))
        left<-anchors(x,type='first');right<-anchors(x,type='second')
        coordinates<-data.frame(chrom1=as.character(seqnames(left)),start1=start(left)-1L,end1=end(left),chrom2=as.character(seqnames(right)),start2=start(right)-1L,end2=end(right))
        write.table(cbind(coordinates,assay(x,'counts')),paste0('counts_',width,'.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
        write.table(cbind(coordinates,offsets),paste0('offsets_',width,'.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
        qc[[as.character(width)]]<-list(resolution=width,bin_pairs=nrow(x),samples=labels,retained_pairs=totals,count_totals=as.numeric(colSums(assay(x,'counts'))),normalization=method,coordinate_system='zero-based half-open')
    }
    write_json(list(method='DIFFHIC_COUNTING_CANDIDATE',historical_equivalence='NOT_CLAIMED',resolutions=qc), 'counting_qc.json',pretty=TRUE,auto_unbox=TRUE)
} else stop('Unknown operation')
writeLines(c(R.version.string,paste('diffHic',packageVersion('diffHic')),paste('csaw',packageVersion('csaw')),paste('edgeR',packageVersion('edgeR'))),'versions.txt')
