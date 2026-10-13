#!/usr/bin/env Rscript
# Render the selected notebook designs; filter changes affect display only.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==9)
settings <- jsonlite::fromJSON(args[8])
code <- readLines(args[1])
ranking <- read.delim(args[2], check.names=FALSE)
stopifnot('factor' %in% names(ranking), nrow(ranking)>=20,
          settings$empty_plot_policy=='center_warning')
dir.create(args[7], recursive=TRUE, showWarnings=FALSE)
quote_r <- function(x) encodeString(x, quote='"')
replace_assignment <- function(script,key,value) {
 idx <- grep(paste0('^',key,'[[:space:]]*<-[[:space:]]*'),script)
 stopifnot(length(idx)==1)
 script[idx] <- paste(key,'<-',value)
 script
}
make_env <- function() {
 env <- new.env(parent=globalenv())
 env$source_label <- args[5]; env$target_label <- args[6]
 env$fc_limits <- settings$limits; env$fc_breaks <- settings$breaks
 env
}
center_panel <- function(tf,out_pdf) {
 # This fallback handles only the selected design's explicit no-edge errors.
 # Recalculate the unfiltered neighbourhood even if filtering failed early.
 suppressPackageStartupMessages({library(ggplot2);library(igraph);library(patchwork)})
 de <- read.delim(args[4],check.names=FALSE)
 row <- de[match(tf,de$gene_id),,drop=FALSE]
 fc <- row$log2FoldChange
 net <- read.delim(args[3],check.names=FALSE)
 graph <- graph_from_data_frame(net[,c('source','target','weight')],directed=TRUE)
 distances_tf <- if(tf %in% V(graph)$name) as.numeric(distances(graph,v=tf,mode='out',weights=NA)) else numeric()
 inf <- ranking[match(tf,ranking$factor),,drop=FALSE]
 warning <- paste0('Warning, ',tf,' edge values fall below the minimum cutoff value of ',settings$min_edge_weight_plot)
 p <- ggplot(data.frame(x=0,y=0,fc=fc),aes(x,y)) +
   geom_point(aes(fill=fc),shape=22,size=13,colour='grey85') +
   geom_label(label=tf,nudge_y=.15,fill='white',label.size=0) +
   scale_fill_gradient2(low='blue3',mid='white',high='red3',midpoint=0,
     limits=settings$limits,breaks=settings$breaks,oob=scales::squish,na.value='grey80',
     name=paste0('log2FC\n(',args[6],' vs ',args[5],')')) +
   coord_cartesian(xlim=c(-1,1),ylim=c(-1,1)) + theme_void() +
   labs(title=paste('2-step branching subnet around',tf),caption=warning) +
   theme(plot.caption=element_text(size=9,hjust=.5),plot.title=element_text(hjust=.5,size=12))
 info <- paste(c(paste0('TF: ',tf),paste0('Direct targets (1-step, total): ',sum(distances_tf==1)),
   paste0('Secondary targets (2-step, total): ',sum(distances_tf==2)),
   paste0('Targets within 2 steps (total): ',sum(distances_tf>=1 & distances_tf<=2)),
   paste0('target_score_scaled: ',signif(inf$target_score_scaled,4)),
   paste0('G_score_scaled: ',signif(inf$G_score_scaled,4)),
   paste0('influence_score (scaled): ',signif(inf$influence_score,4))),collapse='\n')
 panel <- wrap_elements(grid::textGrob(info,x=0,y=.5,just='left',gp=grid::gpar(fontsize=10)))
 ggsave(out_pdf,p|panel,width=15,height=6,device=cairo_pdf)
}
statuses <- list()
empty_errors <- c('No plot edges after min_edge_weight_plot filtering; lower it.',
 'TF not found in filtered plot graph. Lower min_edge_weight_plot.',
 'g_sub0 has 0 edges; lower min_edge_weight_plot.',
 'Outward-edge filter removed all plot edges. Disable keep_outward_edges or lower min_edge_weight_plot.',
 'g_final has 0 plot edges; relax plot filters.')
for(tf in head(ranking$factor,20)) {
 script <- code
 values <- list(influence_file=quote_r(args[2]),diffnet_file=quote_r(args[3]),deseq_file=quote_r(args[4]),
   tf_of_interest=quote_r(tf),out_dir=quote_r(args[7]),
   min_edge_weight_plot=settings$min_edge_weight_plot,padj_cutoff=settings$padj_cutoff,
   abs_lfc_cutoff=settings$abs_lfc_cutoff,
   fc_limits=paste0('c(',paste(settings$limits,collapse=','),')'),
   fc_breaks=paste0('c(',paste(settings$breaks,collapse=','),')'))
 for(key in names(values))script <- replace_assignment(script,key,values[[key]])
 out_pdf <- file.path(args[7],paste0('ANANSE_2step_branch_',tf,'_',args[5],'_to_',args[6],'.pdf'))
 script <- replace_assignment(script,'out_pdf',quote_r(out_pdf))
 writeLines(script,file.path(args[7],paste0(tf,'_executed.R')))
 env <- make_env()
 result <- tryCatch({eval(parse(text=script),envir=env);NULL},error=function(e)e)
 status <- 'SELECTED_FORCE_DIRECTED_PLOT'; reason <- ''
 if(inherits(result,'error')) {
   reason <- conditionMessage(result)
   if(reason %in% empty_errors) {
     center_panel(tf,out_pdf);status <- 'CENTER_ONLY_FILTERED_EDGES'
   } else status <- 'PLOT_ERROR'
 }
 statuses[[length(statuses)+1]] <- data.frame(factor=tf,source=args[5],target=args[6],
   plot_status=status,reason=reason,pdf=basename(out_pdf))
}
write.table(do.call(rbind,statuses),file.path(args[7],'plot_status.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if(any(vapply(statuses,function(x)x$plot_status=='PLOT_ERROR',logical(1))))stop('Unexpected plotting error; inspect plot_status.tsv')
script <- readLines(args[9])
top10_pdf <- file.path(args[7],paste0('ANANSE_top10_',args[5],'_to_',args[6],'.pdf'))
values <- list(influence_file=quote_r(args[2]),diffnet_file=quote_r(args[3]),deseq_file=quote_r(args[4]),
               out_pdf=quote_r(top10_pdf),min_abs_edge_wt=settings$min_edge_weight_plot)
for(key in names(values))script <- replace_assignment(script,key,values[[key]])
writeLines(script,file.path(args[7],'top10_executed.R'))
eval(parse(text=script),envir=make_env())
jsonlite::write_json(settings,file.path(args[7],'plot_settings.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(c(R.version.string,paste('igraph',packageVersion('igraph')),paste('ggraph',packageVersion('ggraph')),
 paste('ggplot2',packageVersion('ggplot2')),paste('Cairo',packageVersion('Cairo'))),'versions.txt')
