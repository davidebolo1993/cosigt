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

#split region
reg<-unlist(strsplit(unlist(strsplit(basename(gtf), ".", fixed=T))[1], "_"))
start<-as.numeric(reg[3])
end<-as.numeric(reg[4])

#process
df<-import(gtf)
genes<-unique(df$gene_name)
count<-0
for (g in genes) {
    subdf<-as.data.frame(df[df$gene_name == g])
    subdf<-head(subdf[order(subdf$end-subdf$start,decreasing=TRUE),],1)[c("seqnames", "start", "end", "gene_name", "gene_type")]
    if (nrow(subdf)>=1) {
        if (subdf$start >= start && subdf$end <= end) {
            count<-count+1
            subdf$seqnames<-path
            fwrite(subdf, file = bed, sep = "\t", row.names = FALSE, col.names = FALSE, append = T)
        }
    }
}

if (count == 0) {
    subdf<-data.frame(seqnames=path, start=start, end=end, gene_name="unknown", gene_type="unknown")
    fwrite(subdf, file = bed, sep = "\t", row.names = FALSE, col.names = FALSE)
}