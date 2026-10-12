suppressPackageStartupMessages({library(diffHic);library(csaw);library(data.table);library(jsonlite)})
root<-'/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow'
source(file.path(root,'source_v23/bin/fixed_grid_helpers.R'))
# Hand-derived positions cross a 100-bp boundary on the negative strand.
stopifnot(identical(fixed_5prime_bins(c(100L,100L,101L),c(11L,-11L,-10L),100L)$bin_start,c(0L,100L,100L)))
folder<-file.path(root,'runs/diffhic_fixed_contacts_fixture_v23/results/diffhic/fixed_grid_contacts')
counts<-fread(file.path(folder,'fixed_grid_counts.tsv'));offsets<-fread(file.path(folder,'fixed_grid_offsets.tsv'))
stopifnot(identical(as.integer(colSums(counts[,.(State_1_HiC_1,State_2_HiC_1)])),c(7L,7L)),all(counts$start1%%1000==0),all(counts$start2%%1000==0),all(is.finite(as.matrix(offsets[,.(State_1_HiC_1,State_2_HiC_1)]))))
contacts<-lapply(c('State_1','State_2'),function(state){con<-gzfile(file.path(folder,'state_contacts',state,'chrtest_contig.txt.gz'),'rt');on.exit(close(con));read.delim(con,header=FALSE)})
stopifnot(identical(contacts[[1]],contacts[[2]]),nrow(contacts[[1]])>0,all(contacts[[1]][,1]%%1000==0),all(contacts[[1]][,2]%%1000==0),all(contacts[[1]][,1]<=contacts[[1]][,2]),all(is.finite(contacts[[1]][,3])),all(contacts[[1]][,3]>0))
provenance<-fromJSON(file.path(folder,'fixed_contacts_provenance.json'))
stopifnot(provenance$method=='diffhic_loess_fixed_grid',all(provenance$fixed_grid_count_totals==7))
out<-file.path(root,'evidence/fixed_contacts_fixture_v01');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_json(list(status='PASS_FIXED_GRID_RETAINED_PAIRS',retained_pair_counts=c(7L,7L),negative_strand_boundary='PASS_HAND_DERIVED_5PRIME_COORDINATES',grid_alignment='PASS',upper_triangle='PASS',state_contacts_identical='PASS_IDENTICAL_SYNTHETIC_LIBRARIES',offsets_finite=TRUE,contact_values_finite_positive=TRUE,method='Synthetic qualification, not historical equivalence'),file.path(out,'qualification.json'),pretty=TRUE,auto_unbox=TRUE)
