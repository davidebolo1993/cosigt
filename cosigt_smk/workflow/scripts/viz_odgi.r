#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(rjson)
library(RColorBrewer)
library(scales)
library(data.table)

setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
coverage_file<-args[1]
node_length_file<-args[2]
cluster_file<-args[3]
#chiara.paleni
medoids<-args[5]

#mod nodes
pad_node_ids <- function(ids) {
  node_nums <- as.integer(sub("node\\.", "", ids))
  sprintf("node.%06d", node_nums)
}

#read inputs
#pangenome coverage over nodes 
coverage_data <- fread(coverage_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coverage_long <- coverage_data %>%
    pivot_longer(cols = -path.name, names_to = "node_id", values_to = "coverage") %>%
    rename(path_name = path.name)
coverage_long$node_id <- pad_node_ids(coverage_long$node_id)

#node lengths
node_length_data <- fread(node_length_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
length_df<-data.frame(
    node_id = node_length_data$V1,
    length = node_length_data$V2,
    stringsAsFactors = FALSE
  )
length_df$node_id<-pad_node_ids(length_df$node_id)

#clustering info
clustering_data<-fromJSON(file = cluster_file)
clustering_df <- data.frame(
    path_name = names(clustering_data),
    cluster = as.character(clustering_data),
    stringsAsFactors = FALSE
)

#merge
viz_data <- merge(coverage_long, length_df, by = "node_id")
viz_data <- merge(viz_data, clustering_df, by = "path_name", all.x = TRUE)

viz_data <- viz_data %>%
  arrange(node_id, path_name) %>% 
  group_by(path_name) %>%
  mutate(
    cumulative_length = cumsum(length),
    start_pos = lag(cumulative_length, default = 0),
    end_pos = cumulative_length
  ) %>%
  ungroup()


spectral_colors <- brewer.pal(11, "Spectral")  # 11-color spectral
max_cov_palette <- 2 + length(spectral_colors) - 1  # max coverage before clamping

color_map <- c(
  "0" = "white",
  "1" = "grey60"
)

for (i in 2:max_cov_palette) {
  color_map[as.character(i)] <- spectral_colors[i - 1]  # offset by 1
}

#mod viz data
viz_data$coverage_clamped <- as.character(
  ifelse(viz_data$coverage >= max_cov_palette, max_cov_palette, viz_data$coverage)
)


viz_data$cluster_num<-as.numeric(gsub("HaploGroup","",viz_data$cluster))
clusters<-unique(viz_data[,c("cluster_num","cluster")])
viz_data$cluster<-factor(viz_data$cluster,levels=clusters$cluster)
p <- ggplot(viz_data, aes(x = start_pos, xend = end_pos, 
                          y = path_name, 
                          color = coverage_clamped)) +
  geom_segment(linewidth = 5) +
  #Andrea asked for this
  scale_color_manual(
    values = color_map,
    na.value = "transparent"
  ) +
  scale_y_discrete(
    expand = expansion(add=c(1,1))
  ) +
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Genomic Position",
    y = "Haplotypes",
    color = "Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position="none"
  )

ggsave(args[4], width=30, height=max(5, 0.25*length(unique(viz_data$path_name))), limitsize=FALSE)

#added by chiara.paleni - viz of representative haplotypes
medoids<-fread(args[5],header = F)
viz_medoid<-viz_data[viz_data$path_name %in% medoids$V2,]
unique_mol<-viz_data[!duplicated(viz_data$path_name),]
abundance<-data.frame(table(unique_mol$cluster))
viz_medoid<-merge(viz_medoid,abundance,by.x="cluster",by.y="Var1")
viz_medoid$label<-paste0("HaploGroup",viz_medoid$cluster_num,"\n",viz_medoid$Freq," haplotype(s)")
viz_medoid<-viz_medoid[order(viz_medoid$Freq,decreasing = T),]
viz_medoid$label<-factor(viz_medoid$label,levels=unique(viz_medoid$label))
p <- ggplot(viz_medoid, 
            aes(x = start_pos, xend = end_pos, 
                          y = path_name, 
                          color = coverage_clamped)) +
  geom_segment(linewidth = 5) +
  #Andrea asked for this
  scale_color_manual(
    values = color_map,
    na.value = "transparent"
  ) +
  scale_y_discrete(
    expand = expansion(add=c(1,1))
  ) +
  facet_grid(label ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Genomic Position",
    y = "Haplotypes",
    color = "Coverage"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position="none"
  )

plot_height <- max(5, 0.4*length(unique(medoids$V1)))
ggsave(gsub(".png",".representative.png",args[4]), height=plot_height, width=25, limitsize = FALSE)
