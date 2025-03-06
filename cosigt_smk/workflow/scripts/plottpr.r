#!/usr/bin/Rscript

library(ggplot2)
library(rjson)
library(data.table)
library(scales)

# Main script logic
args <- commandArgs(trailingOnly = TRUE)

# Set data.table threads
setDTthreads(1)

#Helper function #1
get_region <- function(x) {
  tail(unlist(strsplit(x, "/")), 2)[1]
}

#Helper function #2
process_file <- function(file_path, clusters, distances) {
  df <- fread(file_path, header = TRUE)
  sample <- basename(dirname(dirname(file_path)))
  sample<-unlist(strsplit(sample, ".", fixed=TRUE))[1]
  haps <- names(clusters)[grepl(sample, names(clusters))]
  if (length(haps) != 2) {
    true_clusters<-c("ambiguous", "ambiguous")
  } else {
    true.cl.1<-ifelse(clusters[[haps[1]]] == "*", "unclustered", clusters[[haps[1]]])
    true.cl.2<-ifelse(clusters[[haps[2]]] == "*", "unclustered", clusters[[haps[2]]])
    true_clusters<-sort(c(true.cl.1, true.cl.2))
  }
  pred.cl.1<-ifelse(df$cluster.1[1] == "*", "unclustered", df$cluster.1[1])
  pred.cl.2<-ifelse(df$cluster.2[1] == "*", "unclustered", df$cluster.2[1])
  predicted_clusters<-sort(c(pred.cl.1, pred.cl.2))
  predicted_haplos <- sort(c(df[[1]][1], df[[2]][1]))

  list(
    #for table/plotting
    lenient = identical(true_clusters, predicted_clusters),
    strict = length(grep(sample, predicted_haplos)) == 2,
    w_lenient = ifelse(identical(true_clusters, predicted_clusters), "", sample),
    w_strict = ifelse(length(grep(sample, predicted_haplos)) == 2, "", sample),
    g_lenient = ifelse(identical(true_clusters, predicted_clusters), sample, ""),
    g_strict = ifelse(length(grep(sample, predicted_haplos)) == 2, sample, ""),
    #for table
    sample.id = sample,
    hap.1.pred = predicted_haplos[1],
    hap.2.pred = predicted_haplos[2],
    cl.1.pred = predicted_clusters[1],
    cl.2.pred = predicted_clusters[2],
    cl.1.true = true_clusters[1],
    cl.2.true = true_clusters[2],
    cl.1.pred.true.dist = ifelse(true_clusters[1] == "ambiguous" || true_clusters[1] == "unclustered" || predicted_clusters[1] == "unclustered", "-", distances[which(distances$h.group == true_clusters[1]),][[predicted_clusters[1]]]),
    cl.2.pred.true.dist = ifelse(true_clusters[2] == "ambiguous" || true_clusters[2] == "unclustered" || predicted_clusters[2] == "unclustered", "-", distances[which(distances$h.group == true_clusters[2]),][[predicted_clusters[2]]])
  )
}


files <- sort(list.files(args[1], pattern = "sorted_combos.tsv", full.names = TRUE, recursive = TRUE))
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

  btab<-data.frame(sample=sample, hap.1.pred = hap.1.pred, hap.2.pred = hap.2.pred, cl.1.pred = cl.1.pred, cl.2.pred=cl.2.pred, cl.1.true=cl.1.true, cl.2.true=cl.2.true, cl.1.pred.true.dist = cl.1.pred.true.dist, cl.2.pred.true.dist=cl.2.pred.true.dist)
  fwrite(btab, gsub("pdf", paste0(r,".wdist.tsv"), args[3]),col.names=T, row.names=F, sep="\t")

  data.frame(
    category = rep(c("lenient", "strict"), each = 2),
    measure = rep(c("tp", "fn"), 2),
    region = r,
    score = c(sum(lenient_scores), length(region_files) - sum(lenient_scores),
              sum(strict_scores), length(region_files) - sum(strict_scores)),
    ids=c(paste(lenient_g_names[lenient_g_names!=""],collapse=";"),paste(lenient_w_names[lenient_w_names!=""],collapse=";"),paste(strict_g_names[strict_g_names!=""],collapse=";"),paste(strict_w_names[strict_w_names!=""],collapse=";"))
  )

})

tpr_df <- rbindlist(tpr_list)
tpr_df$measure <- factor(tpr_df$measure, levels = c("fn", "tp"))
tpr_df$category <- factor(tpr_df$category, levels = c("strict", "lenient"))

p <- ggplot(tpr_df, aes(x = category, y = score, fill = measure)) +
  geom_bar(position = "fill", stat = "identity", width = 0.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("tp" = "darkred", "fn" = "darkblue")) +
  facet_wrap(~region, ncol = 6) +
  ylab("% score") +
  xlab("region")

ggsave(args[3], width = 30, height=40)
fwrite(tpr_df, gsub("pdf", "tsv", args[3]),col.names=T, row.names=F, sep="\t")

tpr_df$region<-factor(tpr_df$region, levels=unique(sort(tpr_df$region)))
p2<-ggplot(tpr_df, aes(x = region, y = score, fill = measure)) +
  geom_bar(position = "fill", stat = "identity", width = 0.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("tp" = "darkred", "fn" = "darkblue")) +
  facet_wrap(~category, nrow = 2)

ggsave(gsub("pdf", "mod.pdf", args[3]), width = 40, height=20)
