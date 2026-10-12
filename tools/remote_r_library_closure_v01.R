.libPaths(c('/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1',.libPaths()))
roots<-c('data.table','SummarizedExperiment','csaw','readr','dplyr','tibble','stringr','igraph','ggraph','ggplot2','scales','units','Cairo','patchwork')
ip<-installed.packages()
deps<-tools::package_dependencies(roots,db=ip,which=c('Depends','Imports','LinkingTo'),recursive=TRUE)
selected<-unique(c(roots,unlist(deps)))
selected<-selected[selected %in% rownames(ip)]
selected<-selected[is.na(ip[selected,'Priority']) | ip[selected,'Priority']!='base']
write.table(ip[selected,c('Package','Version','LibPath','Built')],row.names=FALSE,quote=FALSE,sep='\t')
