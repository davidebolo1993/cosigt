#!/usr/bin/Rscript

library(ggplot2)
library(rjson)
library(data.table)
library(scales)
library(dplyr)
library(tidyr)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

get_region <- function(x) {
  tail(unlist(strsplit(x, "/")), 2)[1]
}

process_file <- function(file_path, clusters, distances) {

  df <- fread(file_path, header = TRUE)
  sample <- basename(dirname(dirname(file_path)))
  #in case we have .final, then
  sample<-unlist(strsplit(sample, ".", fixed=TRUE))[1]
  haps <- names(clusters)[grepl(sample, names(clusters))]
  if (length(haps) != 2) {
    true_clusters<-c("ambiguous", "ambiguous")
  } else {
    true.cl.1<-clusters[[haps[1]]]
    true.cl.2<-clusters[[haps[2]]]
    true_clusters<-sort(c(true.cl.1, true.cl.2))
  }
  pred.cl.1<- df$cluster.1[1]
  pred.cl.2<- df$cluster.2[1]
  predicted_clusters<-sort(c(pred.cl.1, pred.cl.2))
  predicted_haplos <- sort(c(df[[1]][1], df[[2]][1]))

  list(
    lenient = identical(true_clusters, predicted_clusters),
    strict = length(grep(sample, predicted_haplos)) == 2,
    w_lenient = ifelse(identical(true_clusters, predicted_clusters), "", sample),
    w_strict = ifelse(length(grep(sample, predicted_haplos)) == 2, "", sample),
    g_lenient = ifelse(identical(true_clusters, predicted_clusters), sample, ""),
    g_strict = ifelse(length(grep(sample, predicted_haplos)) == 2, sample, ""),
    sample.id = sample,
    hap.1.pred = predicted_haplos[1],
    hap.2.pred = predicted_haplos[2],
    cl.1.pred = predicted_clusters[1],
    cl.2.pred = predicted_clusters[2],
    cl.1.true = true_clusters[1],
    cl.2.true = true_clusters[2],
    cl.1.pred.true.dist = ifelse(true_clusters[1] == "ambiguous", "-", distances[which(distances$h.group == true_clusters[1]),][[predicted_clusters[1]]]),
    cl.2.pred.true.dist = ifelse(true_clusters[2] == "ambiguous", "-", distances[which(distances$h.group == true_clusters[2]),][[predicted_clusters[2]]])
  )

}

#only consider files with masks
files <- sort(list.files(args[1], pattern = "sorted_combos.tsv", full.names = TRUE, recursive = TRUE))
files<-files[grep("mask", files, invert=T)]
regions <- unique(sapply(files, get_region))

tpr_list <- lapply(regions, function(r) {

    json_file <- file.path(args[2], paste0(r, ".clusters.json"))
    clusters <- fromJSON(file = json_file)
    dist_file <-file.path(args[2], paste0(r, ".clusters.hapdist.tsv"))
    distances<-fread(dist_file)
    region_files <- files[grepl(r, files)]
    results <- lapply(region_files, process_file, clusters = clusters, distances=distances)
    lenient_scores <- sapply(results, `[[`, "lenient")
    strict_scores <- sapply(results, `[[`, "strict")
    lenient_w_names<-sapply(results, `[[`, "w_lenient")
    strict_w_names<-sapply(results, `[[`, "w_strict")
    lenient_g_names<-sapply(results, `[[`, "g_lenient")
    strict_g_names<-sapply(results, `[[`, "g_strict")
    sample<-sapply(results, `[[`, "sample.id")
    hap.1.pred<-sapply(results, `[[`, "hap.1.pred")
    hap.2.pred<-sapply(results, `[[`, "hap.2.pred")
    cl.1.pred<-sapply(results, `[[`, "cl.1.pred")
    cl.2.pred<-sapply(results, `[[`, "cl.2.pred")
    cl.1.true<-sapply(results, `[[`, "cl.1.true")
    cl.2.true<-sapply(results, `[[`, "cl.2.true")
    cl.1.pred.true.dist<-sapply(results, `[[`, "cl.1.pred.true.dist")
    cl.2.pred.true.dist<-sapply(results, `[[`, "cl.2.pred.true.dist")

    #dissimilarity-table, per region
    diff_table<-fread(file.path(args[3], paste0(r, ".tsv")))
    
    qvlisth1<-rep(0,length(sample))
    qvlisth2<-rep(0,length(sample))
    TPR<-rep("FN",length(sample))

    for (i in 1:length(sample)) {
        TPR[i]<-ifelse(cl.1.pred.true.dist[i] == 0 && cl.2.pred.true.dist[i] == 0, "TP", "FN")
        sample_id<-sample[i]
        hapst<-unique(grep(sample_id, diff_table$group.a,value=T))

        if (length(hapst) != 2) {
            qvlisth1[i]<- -9999
            qvlisth2[i]<- -9999
            message("Missing results for sample ", sample_id)
            next
        }
        hap1t<-hapst[1]
        hap2t<-hapst[2]
        hap1p<-hap.1.pred[i]
        hap2p<-hap.2.pred[i]
        #combo1 hap1t-hap1p
        h1th1p<-diff_table[(diff_table$group.a == hap1t & diff_table$group.b == hap1p)][['estimated.difference.rate']]
        h2th2p<-diff_table[(diff_table$group.a == hap2t & diff_table$group.b == hap2p)][['estimated.difference.rate']]
        e1<-h1th1p+h2th2p
        #other combo
        h1th2p<-diff_table[(diff_table$group.a == hap1t & diff_table$group.b == hap2p)][['estimated.difference.rate']]
        h2th1p<-diff_table[(diff_table$group.a == hap2t & diff_table$group.b == hap1p)][['estimated.difference.rate']]
        e2<-h1th2p+h2th1p
        if (e1 <= e2) {
            qvlisth1[i]<-h1th1p
            qvlisth2[i]<-h2th2p
        } else {
            qvlisth1[i]<-h1th2p
            qvlisth2[i]<-h2th1p
        }
    }

    btab<-data.frame(sample=sample, region=r, hap.1.pred = hap.1.pred, hap.2.pred = hap.2.pred, cl.1.pred = cl.1.pred, cl.2.pred=cl.2.pred, cl.1.true=cl.1.true, cl.2.true=cl.2.true, cl.1.pred.true.dist = cl.1.pred.true.dist, cl.2.pred.true.dist=cl.2.pred.true.dist, qv.1=qvlisth1, qv.2=qvlisth2, TPR=TPR)
    fwrite(btab, gsub("pdf", paste0(r,".wdist.tsv"), args[4]),col.names=T, row.names=F, sep="\t")
    btab
})

tpr_df <- rbindlist(tpr_list)
tpr_df <- subset(tpr_df, (qv.1 != -9999 & qv.2 != -9999))

data_long <- tpr_df %>%
  select(qv.1, qv.2, TPR, region) %>%
  pivot_longer(cols = c(qv.1, qv.2), names_to = "qv_type", values_to = "qv_value")

#tpr pctg per mask per region
tpr_summary <- data_long %>%
  group_by(region) %>%
  summarise(
    TP_count = sum(TPR == "TP"),
    total_count = n(),
    TPR_pct = TP_count / total_count * 100, 
    .groups = "drop"
  )

data_long <- left_join(data_long, tpr_summary, by = c("region"))

tpr_summary <- data_long %>% 
  group_by(region) %>% 
  summarise(
    TP_count = unique(TP_count),
    total_count = unique(total_count),
    TPR_pct = unique(TPR_pct),
    max_qv_value = max(qv_value)
  ) %>% 
  ungroup()

ggplot(data_long, aes(x = region, y = qv_value)) +
  geom_violin() + 
  geom_jitter(aes(color = TPR), width = 0.4, alpha = 0.7) + 
  geom_text(
    data = tpr_summary,
    aes(label = sprintf("%.1f%% (%d/%d)", TPR_pct, TP_count, total_count), y = max_qv_value), 
    vjust = -0.5, 
    hjust = 0.5,
    size = 3
  ) +
  theme_bw() +
  labs(x = "region", y = "estimated.difference.rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))

ggsave(args[4],width = 3*length(unique(tpr_summary$region)), height = 5, limitsize=FALSE)