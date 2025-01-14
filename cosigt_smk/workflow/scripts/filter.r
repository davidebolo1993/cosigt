#!/usr/bin/env Rscript
library(data.table)
#set data.table threads
setDTthreads(1)

#parse input
args <- commandArgs(trailingOnly = TRUE)
cov<-args[1]
len<-args[2]
flt<-args[3]
mask<-args[4]

#process
dfcov<-fread(cov, header=T)
dflen<-fread(len, header=F)

keep<-data.frame(V1=dflen$V1, V2 = 1)
m<-unlist(strsplit(flt, ":"))
for (i in 2:ncol(dfcov)) {
   id<-paste0("node.", i-1)
   if ((length(m) == 2 && m[1] == "common_filter") || (length(m) == 1 && m[1] == "common_filter")) {
      if (all(diff(dfcov[[id]]) == 0)) {
            keep$V2[i]<-0        
      }
   }
   if ((length(m) == 2 && !is.na(as.numeric(m[2]))) || (length(m) == 1 && !is.na(as.numeric(m[1])))) {
      threshold<-ifelse(length(m) == 2, as.numeric(m[2]), as.numeric(m[1]))
      if (dflen$V2[i-1] <= threshold) {
            keep$V2[i-1]<-0  
      }
   }
}

#store mask
keep$V1 <- NULL
fwrite(keep, file = mask, sep = "\t", row.names = FALSE, col.names = FALSE)