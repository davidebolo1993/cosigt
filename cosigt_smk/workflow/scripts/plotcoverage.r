#!/usr/bin/Rscript

library(rtracklayer)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

bw <- import.bw(args[1])
df<- data.frame(bw)
df$seqnames<-do.call(c,lapply(as.character(df$seqnames), function(x) unlist(strsplit(x,":", fixed=T))[1]))
n<-length(unique(df$seqnames))
p<-ggplot(df, aes(x=start, y=score, group=seqnames)) +
    geom_point(show.legend=F) + 
    geom_smooth(show.legend=F) + 
    facet_wrap(~seqnames, ncol=10, scales="free_x") + #up to 10 columns - no more 
    theme_linedraw()
ggsave(args[2],width=30, height=3*round(n/10))