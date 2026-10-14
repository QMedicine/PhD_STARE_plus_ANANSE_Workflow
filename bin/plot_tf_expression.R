#!/usr/bin/env Rscript
# A fixed 90-mm vector panel; all labels and shared colour scales are inputs.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
tab <- read.delim(args[1],check.names=FALSE,na.strings=c('NA',''))
cfg <- jsonlite::fromJSON(args[2])
stopifnot(nrow(tab)==20,identical(as.integer(tab$rank),1:20),cfg$width_mm==90)
if (!is.null(cfg$font_file)) {
  stopifnot(file.exists(cfg$font_file))
  font_dir <- file.path(getwd(),'arial_fonts')
  dir.create(font_dir,showWarnings=FALSE)
  file.copy(cfg$font_file,file.path(font_dir,basename(cfg$font_file)),overwrite=TRUE)
  config <- paste0('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig><include ignore_missing="yes">/etc/fonts/fonts.conf</include><dir>',normalizePath(font_dir),'</dir><cachedir>',getwd(),'/font_cache</cachedir></fontconfig>')
  writeLines(config,'fontconfig.xml')
  Sys.setenv(FONTCONFIG_FILE=normalizePath('fontconfig.xml'))
}
matched <- system2('fc-match',c('-f',shQuote('%{family}'),shQuote('Arial')),stdout=TRUE)
if (!any(grepl('Arial',matched,fixed=TRUE))) stop('Arial font is required; set heatmap_arial_font to an accessible Arial font file')
library(grid)
grDevices::cairo_pdf(args[3],width=90/25.4,height=118/25.4,family='Arial',pointsize=6)
grid.newpage()
# All coordinates are measured in mm; labels and both legends remain inside 90 mm.
txt <- function(label,x,y,size=6,just='centre',rot=0) {
  grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,rot=rot,
            gp=gpar(fontfamily='Arial',fontsize=size))
}
txt('STARE + ANANSE Results',45,114,7)
txt('TF Rank and RNA expression',45,110,6)
left <- 29; top <- 103; row_h <- 3.6; cell_w <- 10
expression_palette <- grDevices::colorRamp(c('#ffffff','#005b36'))
fc_palette <- grDevices::colorRamp(c('#2166ac','#ffffff','#e31a1c'))
colour <- function(value,limits,palette) {
  if(!is.finite(value))return('#d9d9d9')
  if(identical(palette,fc_palette)) {
    scaled <- if(value<0) .5+.5*max(value,limits[1])/abs(limits[1]) else .5+.5*min(value,limits[2])/limits[2]
  } else scaled <- (min(max(value,limits[1]),limits[2])-limits[1])/diff(limits)
  rgb(palette(scaled)/255)
}
for(i in 1:20) {
  y <- top-(i-.5)*row_h
  txt(paste0(i,'. ',tab$factor[i]),left-1,y,6,just='right')
  vals <- c(tab$source_log2_TPM_plus1[i],tab$target_log2_TPM_plus1[i],tab$log2FoldChange[i])
  for(j in 1:3) {
    fill <- if(j<3)colour(vals[j],cfg$tpm_limits,expression_palette) else colour(vals[j],cfg$fc_limits,fc_palette)
    grid.rect(x=unit(left+(j-.5)*cell_w,'mm'),y=unit(y,'mm'),width=unit(cell_w,'mm'),height=unit(row_h,'mm'),gp=gpar(fill=fill,col='#555555',lwd=.25))
  }
}
labels <- c(cfg$source,cfg$target,'Target vs source logFC')
for(j in 1:3)txt(labels[j],left+(j-.5)*cell_w,top-20*row_h-2,5.5,just='right',rot=50)
legend <- function(title,x,y,limits,palette,diverging=FALSE) {
  txt(title,x,y+24,5.5,just='left')
  n <- 100
  for(i in 1:n) {
    value <- limits[1]+(i-.5)/n*diff(limits)
    grid.rect(x=unit(x+1.5,'mm'),y=unit(y+(i-.5)/n*20,'mm'),width=unit(3,'mm'),height=unit(20/n+.02,'mm'),gp=gpar(fill=colour(value,limits,palette),col=NA))
  }
  for(value in unique(c(limits[1],if(diverging)0,limits[2])))txt(format(value,trim=TRUE),x+5,y+(value-limits[1])/diff(limits)*20,5.5,just='left')
}
legend('log2(TPM+1)',65,73,cfg$tpm_limits,expression_palette)
legend('log2FC',65,39,cfg$fc_limits,fc_palette,TRUE)
txt('Grey: missing RNA value',65,30,5,just='left')
txt('Rows follow influence rank; colour limits shared across contrasts',45,3,5)
dev.off()
