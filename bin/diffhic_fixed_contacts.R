#!/usr/bin/env Rscript
# Fixed-grid contacts from retained restriction-filtered read pairs, without bin projection.
suppressPackageStartupMessages({library(diffHic);library(csaw);library(data.table);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=7)
param<-readRDS(args[1]);sheet<-read.delim(args[2],stringsAsFactors=FALSE,check.names=FALSE)
width<-as.integer(args[3]);restrict<-strsplit(args[4],',',fixed=TRUE)[[1]];restrict<-restrict[nzchar(restrict)]
method<-args[5];source(args[6]);files<-args[-seq_len(6)]
if(width<=0 || method!='loess_offsets')stop('Unsupported fixed-grid contact settings')
hic<-sheet[sheet$assay=='HiC',]
hic$group<-if('pool_group' %in% names(hic))ifelse(!is.na(hic$pool_group) & nzchar(hic$pool_group),hic$pool_group,hic$sample_name) else hic$sample_name
groups<-unique(hic$group);labels<-sub('\\.trimmed\\.h5$','',basename(files));order<-match(groups,labels)
if(anyNA(order)||length(order)!=length(files)||anyDuplicated(labels))stop('Fixed-grid library/index identities differ')
files<-files[order];labels<-labels[order]
states<-vapply(labels,function(label){s<-unique(hic$cell_state[hic$group==label]);if(length(s)!=1)stop('Conflicting pooled state');s},character(1))
if(length(files)<2)stop('Joint offsets require two libraries')
lengths<-seqlengths(param$fragments)
if(anyNA(lengths))stop('Reference contig lengths are missing')
if(length(restrict) && !all(restrict %in% names(lengths)))stop('Requested contact chromosomes are absent from reference')
chromosomes<-if(length(restrict))restrict else names(lengths)
pieces<-list();retained<-integer(length(files));total_before_restrict<-integer(length(files));piece_index<-1L
for(lib in seq_along(files)){
    combinations<-loadChromos(files[lib])
    for(i in seq_len(nrow(combinations))){
        c1<-as.character(combinations$anchor1[i]);c2<-as.character(combinations$anchor2[i])
        reads<-loadData(files[lib],c1,c2);total_before_restrict[lib]<-total_before_restrict[lib]+nrow(reads)
        if(!(c1 %in% chromosomes && c2 %in% chromosomes))next
        # Positions are 1-based leftmost coordinates, with signed alignment lengths.
        # Negative-strand 5-prime endpoints are the rightmost aligned base.
        first<-fixed_5prime_bins(reads$anchor1.pos,reads$anchor1.len,width)
        second<-fixed_5prime_bins(reads$anchor2.pos,reads$anchor2.len,width)
        p1<-first$endpoint0;p2<-second$endpoint0
        if(any(p1<0L|p2<0L|p1>=lengths[c1]|p2>=lengths[c2]))stop('Retained endpoint is outside reference')
        b1<-first$bin_start;b2<-second$bin_start
        if(c1==c2){low<-pmin(b1,b2);b2<-pmax(b1,b2);b1<-low}
        if(match(c1,chromosomes)>match(c2,chromosomes)){
            tmp<-c1;c1<-c2;c2<-tmp;tmp<-b1;b1<-b2;b2<-tmp
        }
        tab<-data.table(chrom1=c1,start1=b1,chrom2=c2,start2=b2)[,.(count=.N),by=.(chrom1,start1,chrom2,start2)]
        tab[,sample:=labels[lib]];pieces[[piece_index]]<-tab;piece_index<-piece_index+1L;retained[lib]<-retained[lib]+nrow(reads)
    }
}
if(any(retained==0L)||!length(pieces))stop('No retained fixed-grid interactions for a library')
long<-rbindlist(pieces);wide<-dcast(long,chrom1+start1+chrom2+start2~sample,value.var='count',fill=0L)
setorder(wide,chrom1,start1,chrom2,start2)
mat<-as.matrix(wide[,..labels]);storage.mode(mat)<-'integer'
if(!identical(as.numeric(colSums(mat)),as.numeric(retained)))stop('Fixed-grid count conservation failed')
se<-SummarizedExperiment(assays=list(counts=mat));colData(se)$totals<-retained
se<-csaw::normOffsets(se)
cpm<-calculateCPM(se,use.offsets=TRUE,use.norm.factors=FALSE,log=FALSE,assay.id='counts')
if(any(!is.finite(cpm))||any(cpm<0)||any(!is.finite(assay(se,'offset'))))stop('Invalid fixed-grid normalized contacts')
fwrite(wide,'fixed_grid_counts.tsv',sep='\t')
fwrite(cbind(wide[,.(chrom1,start1,chrom2,start2)],as.data.table(assay(se,'offset'))),'fixed_grid_offsets.tsv',sep='\t')
records<-list()
for(state in unique(states)){
    columns<-which(states==state);values<-rowMeans(cpm[,columns,drop=FALSE]);folder<-file.path('state_contacts',state);dir.create(folder,recursive=TRUE)
    contacts<-list()
    for(chrom in chromosomes){
        rows<-which(wide$chrom1==chrom & wide$chrom2==chrom & values>0)
        if(!length(rows)){if(length(restrict))stop('No positive cis contacts for required chromosome ',chrom,' in ',state);next}
        path<-file.path(folder,paste0('chr',sub('^chr','',chrom),'.txt.gz'))
        if(any(wide$start1[rows]%%width!=0L|wide$start2[rows]%%width!=0L))stop('Non-grid contact start')
        con<-gzfile(path,'wt');write.table(data.frame(start1=wide$start1[rows],start2=wide$start2[rows],value=values[rows]),con,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE);close(con)
        contacts[[chrom]]<-length(rows)
    }
    records[[state]]<-list(cell_state=state,samples=labels[columns],cis_contacts=contacts,summary='arithmetic mean of normalized library CPM',binning='retained 5-prime endpoint, zero-based fixed grid',resolution=width)
}
write_json(list(method='diffhic_loess_fixed_grid',historical_equivalence='NOT_CLAIMED',normalization=method,chosen_res=width,library_order=labels,retained_pairs_before_chromosome_restriction=total_before_restrict,retained_pairs_after_chromosome_restriction=retained,fixed_grid_count_totals=as.numeric(colSums(mat)),coordinates='zero-based bin starts; canonical upper triangle; no rounded-bin projection',states=records),'fixed_contacts_provenance.json',pretty=TRUE,auto_unbox=TRUE)
writeLines(c(R.version.string,paste('diffHic',packageVersion('diffHic')),paste('csaw',packageVersion('csaw')),paste('data.table',packageVersion('data.table'))),'versions.txt')
