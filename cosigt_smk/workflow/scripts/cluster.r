#!/usr/bin/env Rscript 

#libraries
library(data.table)
library(reshape2)
library(NbClust)
library(rjson)
library(dendextend)
library(ggplot2)

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
#we can add min number of clusters and max number of clusters here
res<-NbClust(diss=distanceMatrix, method = "average", index ="silhouette", distance=NULL,max.nc=round(length(unique(df$group.a))/3))$Best.partition
res.list <- lapply(split(res, names(res)), unname)
named_res<-lapply(res.list, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
jout<-toJSON(named_res)
#write json out
write(jout, args[2])
#reverse json
reversed_data <- list()
for (key in names(jout)) {
  value <- jout[[key]]
  if (!is.null(reversed_data[[value]])) {
    reversed_data[[value]] <- c(reversed_data[[value]], key)
  } else {
    reversed_data[[value]] <- key
  }
}
haplotable <- data.frame(
  haplotype.name = unlist(reversed_data),
  haplotype.group = rep(names(reversed_data), lengths(reversed_data))
)
rownames(haplotable)<-c()
fwrite(haplotable, filename=gsub(".json", ".tsv", args[2]), row.names=F,col.names=T,sep="\t")

#simplify names, in distMatrix for hclust
colnames(regularMatrix)<-do.call(c,lapply(colnames(regularMatrix), function(x) (unlist(strsplit(x,":"))[1])))
rownames(regularMatrix)<-do.call(c,lapply(rownames(regularMatrix), function(x) (unlist(strsplit(x,":"))[1])))
distanceMatrix <- as.dist(regularMatrix)
k<-as.numeric(tail(sort(res),1))

#plot
hc <- hclust(distanceMatrix, method = "average")
dend <- as.dendrogram(hc) %>%set("branches_k_color" ,k=k) %>% set("labels_cex", .6) 
ggd1 <- as.ggdend(dend)
p<-ggplot(ggd1,horiz = TRUE) + theme_minimal()
ggsave(gsub(".json", ".pdf", args[2]), width=20,height=15)
