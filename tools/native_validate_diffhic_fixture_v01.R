suppressPackageStartupMessages({library(diffHic);library(jsonlite)})
root<-'/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow'
folder<-file.path(root,'runs/diffhic_count_fixture_v21/results/diffhic/counts')
records<-list()
fragments<-cutGenome(file.path(root,'inputs/hic_presplit_fixture_v01/reference.fa'),'GATC',4)
for(width in c(500L,1000L)){
    x<-readRDS(file.path(folder,paste0('interactions_',width,'.rds')))
    counts<-assay(x,'counts');offsets<-assay(x,'offset')
    stopifnot(nrow(counts)>0,ncol(counts)==2,all(is.finite(offsets)),all(counts>=0),identical(counts[,1],counts[,2]),identical(colnames(counts),c('State_1_HiC_1','State_2_HiC_1')))
    totals<-as.numeric(colSums(counts));stopifnot(all(totals>0),all(totals==as.numeric(x$totals)))
    table<-read.delim(file.path(folder,paste0('counts_',width,'.tsv')),check.names=FALSE)
    stopifnot(nrow(table)==nrow(counts),identical(as.matrix(table[,colnames(counts)]),counts))
    stopifnot(all(table$start1 %in% (start(fragments)-1L)),all(table$start2 %in% (start(fragments)-1L)),all(table$end1 %in% end(fragments)),all(table$end2 %in% end(fragments)))
    stopifnot(identical(table$start1,as.integer(start(anchors(x,type='first'))-1L)),identical(table$end2,as.integer(end(anchors(x,type='second')))))
    records[[as.character(width)]]<-list(resolution=width,bin_pairs=nrow(counts),retained_pair_counts=totals,offsets_finite=TRUE,sample_order_exact=TRUE,table_values_exact=TRUE,binning='restriction-site-rounded squareCounts; not fixed-grid cooler bins')
}
out<-file.path(root,'evidence/diffhic_count_fixture_v01');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_json(list(status='PASS_SYNTHETIC_TWO_RESOLUTION_COUNTING',resolutions=records,method='diffHic 1.40.0 current-library candidate; synthetic pruning override; no historical or biological equivalence claim'),file.path(out,'qualification.json'),pretty=TRUE,auto_unbox=TRUE)
