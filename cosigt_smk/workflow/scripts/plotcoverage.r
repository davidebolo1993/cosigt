#!/usr/bin/Rscript

library(rtracklayer)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

bw <- import.bw(args[1])
df<- data.frame(bw)
p<-ggplot(df, aes(x=start, y=score, group=seqnames)) +
    geom_point(show.legend=F) + 
    geom_smooth(show.legend=F) + 
    facet_wrap(~seqnames, ncol=5, scales="free_x") +
    theme_linedraw()
ggsave(args[2],width=20, height=20)