#!/usr/bin/Rscript

library(data.table)
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
tpr_file<-args[1]
qv_file<-args[2]
region<-args[3]
df_out<-args[4]
tpr_df<-subset(fread(tpr_file, header=T), tpr != "missing")
qv_df<-fread(qv_file, header=F)
tpr_df$qv.1 <- qv_df$V2[match(tpr_df$sample.id, qv_df$V1)]
tpr_df$qv.2 <- qv_df$V3[match(tpr_df$sample.id, qv_df$V1)]
tpr_df$region<-region
fwrite(tpr_df, df_out, row.names = FALSE, col.names = TRUE, sep = "\t")