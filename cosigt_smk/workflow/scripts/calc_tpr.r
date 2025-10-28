#!/usr/bin/Rscript

library(ggplot2)
library(rjson)
library(data.table)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

#parse input
dissimilarity_table<-args[1]
region<-unlist(strsplit(basename(dissimilarity_table), ".", fixed=T))[1]
clusters_json<-args[2]
clusters_table<-args[3]
output_df<-args[4]
cosigt_predictions<-args[5:length(args)]

#process
process_file <- function(file_path, clusters, cluster_distances, graph_distances) {
    df <- fread(file_path, header = TRUE)
    sample <- basename(dirname(dirname(dirname(file_path)))) #sample name with respect to cosigt new output
    #just in case we have are using .crams with <sample>.final, then
    sample<-unlist(strsplit(sample, ".", fixed=TRUE))[1]
    haps <- sort(names(clusters)[grepl(sample, names(clusters))])
    predicted_haplotypes <- sort(c(df[[1]][1], df[[2]][1]))
    predicted_clusters<-sort(c(df$cluster.1[1], df$cluster.2[1]))
    #only process non-ambiguous for the time being    
    if (length(haps) != 2) {
        list(
            sample.id = sample,
            hap.1.pred = predicted_haplotypes[1],
            hap.2.pred = predicted_haplotypes[2],
            hap.1.true = "missing",
            hap.2.true = "missing",
            cl.1.pred = predicted_clusters[1],
            cl.2.pred = predicted_clusters[2],
            cl.1.true = "missing",
            cl.2.true = "missing",
            cl.1.pred.true.dist= "missing",
            cl.2.pred.true.dist= "missing",
            n.clust=length(unique(unlist(clusters))),
            edr.1="missing",
            edr.2="missing",
            tpr="missing"
        )
    } else {
        true_clusters<-sort(c(clusters[[haps[1]]],clusters[[haps[2]]]))
        #if pred and exp fall in the same cluster, TP, else FN
        cl.1.pred.true.dist = cluster_distances[h.group == true_clusters[1], get(predicted_clusters[1])]
        cl.2.pred.true.dist = cluster_distances[h.group == true_clusters[2], get(predicted_clusters[2])]
        tpr<-ifelse(cl.1.pred.true.dist == 0 && cl.2.pred.true.dist == 0, "TP", "FN")
        h1t<-haps[1]
        h2t<-haps[2]
        h1p<-predicted_haplotypes[1]
        h2p<-predicted_haplotypes[2]
        #h1t = "h1 true"
        #h2t = "h2 true"
        #h1p = "h1 predicted"
        #h2p = "h2 predicted"
        #1st combo
        h1th1p<-graph_distances[group.a == h1t & group.b == h1p, get("estimated.difference.rate")]
        h2th2p<-graph_distances[group.a == h2t & group.b == h2p, get("estimated.difference.rate")]    
        #sum of the 2 errors
        e1<-h1th1p+h2th2p
        #2nd combo
        h1th2p<-graph_distances[group.a == h1t & group.b == h2p, get("estimated.difference.rate")]
        h2th1p<-graph_distances[group.a == h2t & group.b == h1p, get("estimated.difference.rate")]
        e2<-h1th2p+h2th1p
        #haplotyping error - edr, but from graph ~qv?
        if (e1 <= e2) {
            edr.1<-h1th1p
            edr.2<-h2th2p
        } else {
            edr.1<-h1th2p
            edr.2<-h2th1p
        }
        list(
            sample.id = sample,
            hap.1.pred = predicted_haplotypes[1],
            hap.2.pred = predicted_haplotypes[2],
            hap.1.true = haps[1],
            hap.2.true = haps[2],
            cl.1.pred = predicted_clusters[1],
            cl.2.pred = predicted_clusters[2],
            cl.1.true = true_clusters[1],
            cl.2.true = true_clusters[2],
            cl.1.pred.true.dist=cl.1.pred.true.dist,
            cl.2.pred.true.dist=cl.2.pred.true.dist,
            n.clust=length(unique(unlist(clusters))),
            edr.1=edr.1,
            edr.2=edr.2,
            tpr=tpr
        )
    }
}

#load
clusters <- fromJSON(file = clusters_json)
cluster_distances<-fread(clusters_table, header=T)
graph_distances<-fread(dissimilarity_table, header=T)
#process
results <- lapply(cosigt_predictions, process_file, clusters=clusters, cluster_distances=cluster_distances, graph_distances=graph_distances)
#out.table
tpr_df<-rbindlist(results)
#write
fwrite(tpr_df,output_df,col.names=T, row.names=F, sep="\t")