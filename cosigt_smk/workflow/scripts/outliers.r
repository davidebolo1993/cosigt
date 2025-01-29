#!/usr/bin/env Rscript

library(data.table)
library(reshape2)
library(dbscan)

#set data.table threads
setDTthreads(1)

#parse input
args <- commandArgs(trailingOnly = TRUE)
bed_in<-args[1]
bed_out<-args[2]
#read
df<-fread(bed_in)
#process
difflist<-list()
for (i in 1:nrow(df)) {
    h1<-df$V1[i]
    h1len<-df$V3[i]-df$V2[i]
    for (j in 1:nrow(df)) {
        h2<-df$V1[j]
        h2len<-df$V3[j]-df$V2[j]
        pct<-1-min(h1len,h2len)/max(h1len,h2len)
        ddf<-data.frame(h1=h1, h2=h2, diff=pct)
        difflist[[paste0(h1,h2)]]<-ddf
    }
}
diffdf<-do.call(rbind,difflist)
rownames(diffdf)<-NULL
regularMatrix <- acast(diffdf, h1~h2, value.var = "diff")
distanceMatrix <- as.dist(regularMatrix)
res <- dbscan(distanceMatrix, eps=0.2, minPts = 10)$cluster #
#noise points are 0s
outliers<-labels(distanceMatrix)[which(res == 0)]
#remove from df
if (length(outliers) != 0) {
    df<-df[-which(df$V1 %in% outliers)]
}
#write out
fwrite(df, bed_out, row.names = FALSE, col.names = FALSE, sep = "\t")

