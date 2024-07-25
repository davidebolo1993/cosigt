#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(rtracklayer)
library(ggplot2)

files<-sort(list.files(args[1], pattern="*.all.bw", full.names=T,recursive=T))

for (f in files) {

    name<-basename(dirname(f))
    bw <- import.bw(f)
    df<- data.frame(bw)
    p<-ggplot(df, aes(x=start, y=score, group=seqnames)) +
    geom_point(show.legend=F) + 
    geom_smooth(show.legend=F) + 
    facet_wrap(~seqnames, ncol=5, scales="free_x") +
    theme_linedraw()
    ggsave(paste0(f,".pdf"),width=20, height=20)

}