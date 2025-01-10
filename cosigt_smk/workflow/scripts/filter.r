#!/usr/bin/env Rscript
library(data.table)
#set data.table threads
setDTthreads(1)

#parse input
args <- commandArgs(trailingOnly = TRUE)
cov<-args[1]
len<-args[2]
mask<-args[3]

#process
dfcov<-fread(cov, header=T)
dflen<-fread(len, header=F)

keep<-data.frame(V1=dflen$V1, V2 = 1)
for (i in 2:ncol(dfcov)) {
   id<-paste0("node.", i-1)
   if (all(diff(dfcov[[id]]) == 0)) {
        keep$V2[i]<-0        
   }
   if (dflen$V2[i-1] <= 1) {
         keep$V2[i-1]<-0  
   }
}

#store mask
fwrite(keep, file = mask, sep = "\t", row.names = FALSE, col.names = FALSE)