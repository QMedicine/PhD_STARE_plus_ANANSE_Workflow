#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==11,args[11]=='fail')
# args: archived notebook, influence, diffnetwork, H13CORE DE, source, target, output directory
# The generic template is the exact user-selected branching section with supplied labels/paths.
code<-readLines(args[1])
stopifnot(any(grepl('# - Uses FILTERED diffnetwork for PLOT (min_edge_weight)',code,fixed=TRUE)))
ranking<-read.delim(args[2],check.names=FALSE)
stopifnot('factor' %in% names(ranking),nrow(ranking)>=20)
dir.create(args[7],recursive=TRUE)
quote_r<-function(x)encodeString(x,quote='"')
statuses<-list()
for(tf in head(ranking$factor,20)) {
 script<-code
 replacements<-list(influence_file=args[2],diffnet_file=args[3],deseq_file=args[4],tf_of_interest=tf,out_dir=args[7])
 for(key in names(replacements)) {
  idx<-grep(paste0('^',key,'[[:space:]]*<-[[:space:]]*'),script);stopifnot(length(idx)==1)
  script[idx]<-paste(key,'<-',quote_r(replacements[[key]]))
 }
 idx<-grep('^out_pdf <- ',script);stopifnot(length(idx)==1)
 script[idx]<-paste0('out_pdf <- file.path(out_dir, ',quote_r(paste0('ANANSE_2step_branch_',tf,'_',args[5],'_to_',args[6],'.pdf')),')')
 env<-new.env(parent=globalenv());env$source_label<-args[5];env$target_label<-args[6];env$fc_limits<-as.numeric(args[8:9]);env$fc_breaks<-as.numeric(strsplit(args[10],',',fixed=TRUE)[[1]])
 result<-tryCatch({eval(parse(text=script),envir=env);NULL},error=function(e)e)
 status<-'ORIGINAL_FILTERED_PLOT';reason<-''
 if(inherits(result,'error')) {
  reason<-conditionMessage(result)
  # Record the original error without relaxing the notebook filters.
  status<-'ORIGINAL_PLOT_ERROR'
 }
 statuses[[length(statuses)+1]]<-data.frame(factor=tf,source=args[5],target=args[6],plot_status=status,reason=reason,pdf=env$out_pdf)

}
write.table(do.call(rbind,statuses),file.path(args[7],'plot_status.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(R.version.string,paste('igraph',packageVersion('igraph')),paste('ggraph',packageVersion('ggraph')),paste('ggplot2',packageVersion('ggplot2')),paste('Cairo',packageVersion('Cairo'))),'versions.txt')

if(any(vapply(statuses,function(x)x$plot_status=='ORIGINAL_PLOT_ERROR',logical(1))))stop('One or more original branching plots failed; inspect plots/plot_status.tsv. Filters were preserved.')
