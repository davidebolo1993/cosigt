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

# Function to add padding for small clusters and center haplotypes + add facet padding for all clusters
# THIS WAS BUILT WITH THE HELP OF CLAUDE
add_cluster_padding <- function(data, min_size = 1) {
  # Count paths per cluster
  cluster_counts <- data %>%
    group_by(cluster) %>%
    summarise(n_paths = n_distinct(path_name), .groups = "drop")
  
  # Identify small clusters that need centering padding
  small_clusters <- cluster_counts %>%
    filter(n_paths < min_size) %>%
    pull(cluster)
  
  # Get all clusters for facet padding
  all_clusters <- unique(data$cluster)
  
  padding_rows <- list()
  
  # Add facet padding (top/bottom whitespace) for ALL clusters
  for (cluster_name in all_clusters) {
    cluster_data <- filter(data, cluster == cluster_name)
    
    # Get unique nodes for this cluster
    unique_nodes <- unique(cluster_data$node_id)
    node_lengths <- cluster_data %>% 
      select(node_id, length) %>% 
      distinct()
    
    # Add top facet padding
    padding_data_top <- data.frame(
      node_id = unique_nodes,
      path_name = paste0("__FACET_PADDING_TOP__", cluster_name),
      coverage = NA,
      cluster = cluster_name,
      stringsAsFactors = FALSE
    )
    padding_data_top <- merge(padding_data_top, node_lengths, by = "node_id")
    padding_rows[[length(padding_rows) + 1]] <- padding_data_top
    
    # Add bottom facet padding
    padding_data_bottom <- data.frame(
      node_id = unique_nodes,
      path_name = paste0("__FACET_PADDING_BOTTOM__", cluster_name),
      coverage = NA,
      cluster = cluster_name,
      stringsAsFactors = FALSE
    )
    padding_data_bottom <- merge(padding_data_bottom, node_lengths, by = "node_id")
    padding_rows[[length(padding_rows) + 1]] <- padding_data_bottom
  }
  
  # Add centering padding only for small clusters
  for (cluster_name in small_clusters) {
    cluster_data <- filter(data, cluster == cluster_name)
    n_existing <- n_distinct(cluster_data$path_name)
    n_padding_needed <- min_size - n_existing
    
    # Calculate padding distribution to center haplotypes
    padding_before <- floor(n_padding_needed / 2)
    padding_after <- n_padding_needed - padding_before
    
    # Get unique nodes for this cluster
    unique_nodes <- unique(cluster_data$node_id)
    node_lengths <- cluster_data %>% 
      select(node_id, length) %>% 
      distinct()
    
    # Create padding paths BEFORE real haplotypes
    for (i in 1:padding_before) {
      padding_path_name <- paste0("__PADDING_BEFORE__", cluster_name, "__", i)
      
      padding_data <- data.frame(
        node_id = unique_nodes,
        path_name = padding_path_name,
        coverage = NA,
        cluster = cluster_name,
        stringsAsFactors = FALSE
      )
      
      padding_data <- merge(padding_data, node_lengths, by = "node_id")
      padding_rows[[length(padding_rows) + 1]] <- padding_data
    }
    
    # Create padding paths AFTER real haplotypes
    for (i in 1:padding_after) {
      padding_path_name <- paste0("__PADDING_AFTER__", cluster_name, "__", i)
      
      padding_data <- data.frame(
        node_id = unique_nodes,
        path_name = padding_path_name,
        coverage = NA,
        cluster = cluster_name,
        stringsAsFactors = FALSE
      )
      
      padding_data <- merge(padding_data, node_lengths, by = "node_id")
      padding_rows[[length(padding_rows) + 1]] <- padding_data
    }
  }
  
  # Combine original data with padding
  if (length(padding_rows) > 0) {
    padding_df <- do.call(rbind, padding_rows)
    
    # Recalculate cumulative positions for padding data
    padding_df <- padding_df %>%
      arrange(node_id, path_name) %>%
      group_by(path_name) %>%
      mutate(
        cumulative_length = cumsum(length),
        start_pos = lag(cumulative_length, default = 0),
        end_pos = cumulative_length
      ) %>%
      ungroup()
    
    data <- rbind(data, padding_df)
  }
  
  return(data)
}

# Apply padding
viz_data <- add_cluster_padding(viz_data, min_size = 1)

cluster_order <- sort(unique(viz_data$cluster))
viz_data$cluster_num <- match(viz_data$cluster, cluster_order)

cluster_spacing <- 1

path_order_df <- viz_data %>%
  select(cluster, path_name) %>%
  distinct() %>%
  arrange(cluster, path_name) %>%
  group_by(cluster) %>%
  mutate(
    # Custom sorting: facet top padding, centering padding before, real paths, centering padding after, facet bottom padding
    sort_order = case_when(
      grepl("__FACET_PADDING_TOP__", path_name) ~ 0,
      grepl("__PADDING_BEFORE__", path_name) ~ 1,
      grepl("__PADDING_AFTER__", path_name) ~ 3,
      grepl("__FACET_PADDING_BOTTOM__", path_name) ~ 4,
      TRUE ~ 2
    ),
    path_name_clean = case_when(
      grepl("__FACET_PADDING_TOP__", path_name) ~ paste0("_facet_top_", gsub(".*__FACET_PADDING_TOP__(.*)", "\\1", path_name)),
      grepl("__FACET_PADDING_BOTTOM__", path_name) ~ paste0("_facet_bottom_", gsub(".*__FACET_PADDING_BOTTOM__(.*)", "\\1", path_name)),
      grepl("__PADDING_BEFORE__", path_name) ~ paste0("_before_", gsub(".*__PADDING_BEFORE__(.*)__.*", "\\1", path_name)),
      grepl("__PADDING_AFTER__", path_name) ~ paste0("_after_", gsub(".*__PADDING_AFTER__(.*)__.*", "\\1", path_name)),
      TRUE ~ path_name
    )
  ) %>%
  arrange(cluster, sort_order, path_name_clean) %>%
  mutate(path_order_in_cluster = row_number()) %>%
  ungroup()

cluster_sizes <- path_order_df %>%
  group_by(cluster) %>%
  summarise(cluster_size = max(path_order_in_cluster), .groups = "drop") %>%
  arrange(cluster) %>%
  mutate(
    cluster_offset = cumsum(lag(cluster_size + cluster_spacing, default = 0))
  )

path_order_df <- path_order_df %>%
  left_join(cluster_sizes, by = "cluster") %>%
  mutate(y_pos = cluster_offset + path_order_in_cluster)

viz_data <- viz_data %>%
  left_join(path_order_df %>% select(cluster, path_name, y_pos), by = c("cluster", "path_name"))

path_labels <- path_order_df %>%
  filter(!grepl("__PADDING|__FACET_PADDING", path_name)) %>%  # Exclude all padding paths
  select(path_name, y_pos, cluster) %>%
  arrange(y_pos)


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

#actual viz
p <- ggplot(viz_data, aes(x = start_pos, xend = end_pos, 
                           y = y_pos, yend = y_pos, 
                           color = coverage_clamped)) +
    geom_segment(linewidth = 5) +
    #Andrea asked for this
    scale_color_manual(
      values = color_map,
      na.value = "transparent"
    ) +
    scale_y_continuous(
      breaks = path_labels$y_pos,
      labels = path_labels$path_name,
      expand = c(0.02, 0.02)
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

ggsave(args[4], width=30, height=max(5, 0.2*length(unique(viz_data$path_name))), limitsize=FALSE)

#added by chiara.paleni - viz of representative haplotypes
medoids<-fread(args[5],header = F)
viz_medoid<-viz_data[viz_data$path_name %in% medoids$V2,]
unique_mol<-viz_data[!duplicated(viz_data$path_name) & !startsWith(viz_data$path_name,"__FACET_PADDING"),]
abundance<-data.frame(table(unique_mol$cluster))
viz_medoid<-merge(viz_medoid,abundance,by.x="cluster",by.y="Var1")
viz_medoid$label<-paste0("HaploGroup",viz_medoid$cluster_num,"\n",viz_medoid$Freq," haplotype(s)")
viz_medoid<-viz_medoid[order(viz_medoid$Freq,decreasing = T),]
viz_medoid$label<-factor(viz_medoid$label,levels=unique(viz_medoid$label))
p <- ggplot(viz_medoid, 
            aes(x = start_pos, xend = end_pos, 
                          y = y_pos, yend = y_pos, 
                          color = coverage_clamped)) +
  geom_segment(linewidth = 5) +
  #Andrea asked for this
  scale_color_manual(
    values = color_map,
    na.value = "transparent"
  ) +
  scale_y_continuous(
    breaks = path_labels$y_pos,
    labels = path_labels$path_name,
    expand = c(0.02, 0.02)
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
