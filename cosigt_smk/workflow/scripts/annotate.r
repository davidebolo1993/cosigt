#!/usr/bin/env Rscript

library(rtracklayer)
library(data.table)

#set data.table threads
setDTthreads(1)

#parse input
args <- commandArgs(trailingOnly = TRUE)
gtf<-args[1]
path<-args[2]
bed<-args[3]

#process
df<-import(gtf)
genes<-unique(df$gene_name)
for (g in genes) {
    subdf<-as.data.frame(df[df$gene_name == g])
    subdf<-head(subdf[order(subdf$end-subdf$start,decreasing=TRUE),],1)[c("seqnames", "start", "end", "gene_name")]
    if (nrow(subdf)>=1) {
        subdf$seqnames<-path
        fwrite(subdf, file = bed, sep = "\t", row.names = FALSE, col.names = FALSE, append = T)
    }
}