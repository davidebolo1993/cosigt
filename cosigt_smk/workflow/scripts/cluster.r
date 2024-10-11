#!/usr/bin/env Rscript 

#libraries
library(data.table)
library(reshape2)
library(NbClust)
library(rjson)

#args
args <- commandArgs(trailingOnly = TRUE)

#do not exceed
setDTthreads(1)

#execute - do not plot
df<-fread(args[1])
df$jaccard.distance<-1-df$jaccard.similarity
regularMatrix <- acast(df, group.a ~ group.b, value.var = "jaccard.distance")
distanceMatrix <- as.dist(regularMatrix)
#calculate silhouette score and best partition
res<-NbClust(diss=distanceMatrix, method = "average", index ="silhouette", distance=NULL)$Best.partition
named_res<-lapply(res, function(k) paste0("HaploGroup", as.list(res[[k]])))
jout<-toJSON(named_res)
#write json out
write(jout, args[2])
