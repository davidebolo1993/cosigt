#!/usr/bin/Rscript

library(ggplot2)
library(rjson)
library(data.table)
library(scales)

# Main script logic
args <- commandArgs(trailingOnly = TRUE)

#Helper function #1
get_region <- function(x) {
  tail(unlist(strsplit(x, "/")), 2)[1]
}

#Helper function #2
process_file <- function(file_path, clusters) {
  df <- fread(file_path, header = TRUE)
  sample <- basename(dirname(dirname(file_path)))
  haps <- names(clusters)[grepl(sample, names(clusters))]
  true_clusters <- sort(c(clusters[[haps[1]]], clusters[[haps[2]]]))
  predicted_clusters <- sort(c(df$cluster.1[1], df$cluster.2[1]))
  predicted_haplos <- c(df[[1]][1], df[[2]][1])
  
  list(
    lenient = identical(true_clusters, predicted_clusters),
    strict = length(grep(sample, predicted_haplos)) == 2
  )
}


files <- sort(list.files(args[1], pattern = "sorted_combos.tsv", full.names = TRUE, recursive = TRUE))
regions <- unique(sapply(files, get_region))

tpr_list <- lapply(regions, function(r) {
  json_file <- file.path(args[2], paste0(r, ".clusters.json"))
  clusters <- fromJSON(file = json_file)
  region_files <- files[grepl(r, files)]
  
  results <- lapply(region_files, process_file, clusters = clusters)
  lenient_scores <- sapply(results, `[[`, "lenient")
  strict_scores <- sapply(results, `[[`, "strict")
  
  data.frame(
    category = rep(c("lenient", "strict"), each = 2),
    measure = rep(c("tp", "fn"), 2),
    region = r,
    score = c(sum(lenient_scores), length(region_files) - sum(lenient_scores),
              sum(strict_scores), length(region_files) - sum(strict_scores))
  )
})

tpr_df <- rbindlist(tpr_list)
tpr_df$measure <- factor(tpr_df$measure, levels = c("fn", "tp"))
tpr_df$category <- factor(tpr_df$category, levels = c("strict", "lenient"))

p <- ggplot(tpr_df, aes(x = region, y = score, fill = measure)) +
  geom_bar(position = "fill", stat = "identity", width = 0.1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("tp" = "darkred", "fn" = "darkblue")) +
  facet_wrap(~category) +
  ylab("% score") +
  xlab("region")

ggsave(args[3], width = 20)