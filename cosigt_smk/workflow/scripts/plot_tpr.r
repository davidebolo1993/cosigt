#!/usr/bin/Rscript

library(ggplot2)
library(ggrepel)
library(data.table)
library(scales)
library(dplyr)
library(tidyr)
library(cowplot)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

tpr_qv_list<-list()

output_plot_prefix<-args[1]
annot_bed<-args[2]
annot_tsv<-fread(annot_bed)
annot_tsv <- annot_tsv %>%
  mutate(id = paste0(as.character(V1), "_", as.character(V2), "_", as.character(V3)))

tpr_files<-args[3:length(args)]
for (i in 1:length(tpr_files)) {
    tpr_qv_list[[i]]<-fread(tpr_files[i], header=T) #this is already cleaned
}
tpr_df<-do.call(rbind,tpr_qv_list)

#info of interest and convert to long format for ggplot2
data_long <- tpr_df %>%
  select(sample.id ,region, n.clust, edr.1, edr.2, qv.1, qv.2, tpr) %>% 
  pivot_longer(cols = c(edr.1, edr.2, qv.1, qv.2), names_to = "metric", values_to = "metric.values")

data_long$pos<-data_long$region
if (any(annot_tsv$V4 != "unkown")) {
  #match annotation in that case
  data_long$region <- annot_tsv$V4[match(data_long$pos, annot_tsv$id)]
}

#calculate tpr_per_region
tpr_summary <- data_long %>%
  group_by(region,pos) %>%
  summarise(
    tp_count = sum(tpr == "TP"),
    total_count = n(),
    tpr_pct = tp_count / total_count * 100, 
    n_clust = unique(n.clust),
    max_edr_value = max(metric.values[which(metric == "edr.1" | metric == "edr.2")]),
    max_qv_value = max(metric.values[which(metric == "qv.1" | metric == "qv.2")]),
  ) %>% 
  ungroup()

#precompute a label for plotting
tpr_summary$label <- sprintf("tpr=%.1f%% (%d/%d), #cl = %d", 
                             tpr_summary$tpr_pct, 
                             tpr_summary$tp_count/2, #for each sample/haplotype we actually have 2 lines (one w/ )
                             tpr_summary$total_count/2, #same as above
                             tpr_summary$n_clust)

#left join on the original table
data_long <- left_join(data_long, tpr_summary, by = c("region", "pos"))

#max edr and max qv
max_edr <- max(tpr_summary$max_edr_value, na.rm = TRUE)
max_qv <- max(tpr_summary$max_qv_value, na.rm = TRUE)

#flag false-negatives
data_long <- data_long %>%
  group_by(region,pos) %>%
  mutate(
    label_point_edr = ifelse(tpr == "FN", sample.id, ""),
    label_point_qv = ifelse(tpr == "FN", sample.id, "")
  ) %>% ungroup()


num_regions <- length(unique(data_long$region))
#number of columns (max 20 per row)
num_cols <- min(10, num_regions)

#edr and qv - split here
data_long_edr<-subset(data_long, (metric %in% c("edr.1", "edr.2")))
data_long_qv<-subset(data_long, (metric %in% c("qv.1", "qv.2")))

#qv
jitter_width <- 0.15
jitter_height_edr <- 0.01 * max_edr
jitter_height_qv <- 0.01 * max_qv

#plot edr first
data_long_edr <- data_long_edr %>%
  group_by(region, pos, sample.id) %>%
  mutate(
    x_jitter = as.numeric(factor(region)) + runif(1, -jitter_width, jitter_width),
    y_jitter = metric.values + runif(1, -jitter_height_edr, jitter_height_edr)
  ) %>%
  ungroup()

p_edr <- ggplot(data_long_edr, aes(x = region, y = metric.values)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  geom_point( 
    aes(x = x_jitter, y = y_jitter, color = tpr, alpha = tpr),
    size = 3
  ) +
  scale_alpha_manual(
    values = c("FN" = 1.0, "TP" = 0.2),
    guide = "none"
  ) +
  geom_text_repel(
    aes(x = x_jitter, y = y_jitter, label = label_point_edr),
    size = 2.5,
    box.padding = 0.5,
    point.padding = 0.1,
    force = 2,
    segment.size = 0.2,
    max.overlaps = 50
  ) +
  geom_text(
    data = tpr_summary,
    aes(label = label), 
    y = max_edr * 1.03, 
    vjust = -0.5, 
    hjust = 0.5,
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, max_edr * 1.05),
    oob = scales::oob_keep
  )+
  facet_wrap(~ pos, ncol = num_cols, scales = "free_x") +
  theme_bw() +
  labs(y = "estimated.difference.rate", x = "") +
  theme(
    axis.ticks.x = element_blank()
  )

#calculate dimensions
plot_width <- min(3 * num_cols, 30)
plot_height <- 5 * ceiling(num_regions / num_cols)
ggsave(paste0(output_plot_prefix, ".edr.png"), plot = p_edr, width = plot_width, height = plot_height, limitsize = FALSE)

#plot qv then
data_long_qv <- data_long_qv %>%
  group_by(region, pos, sample.id) %>%
  mutate(
    x_jitter = as.numeric(factor(region)) + runif(1, -jitter_width, jitter_width),
    y_jitter = metric.values + runif(1, -jitter_height_qv, jitter_height_qv)
  ) %>%
  ungroup()

p_qv <- ggplot(data_long_qv, aes(x = region, y = metric.values)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  geom_point( 
    aes(x = x_jitter, y = y_jitter, color = tpr, alpha = tpr),
    size = 3
  )  +
  scale_alpha_manual(
    values = c("FN" = 1.0, "TP" = 0.2),
    guide = "none"
  ) +
  geom_text_repel(
    aes(x = x_jitter, y = y_jitter, label = label_point_qv),
    size = 2.5,
    box.padding = 0.5,
    point.padding = 0.1,
    force = 2,
    segment.size = 0.2,
    max.overlaps = 50
  ) +
  geom_text(
    data = tpr_summary,
    aes(label = label), 
    y = max_qv * 1.03, 
    vjust = -0.5, 
    hjust = 0.5,
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, max_qv * 1.05),
    oob = scales::oob_keep
  )+
  facet_wrap(~ pos, ncol = num_cols, scales = "free_x") +
  theme_bw() +
  labs(y = "QV", x = "") +
  theme(
    axis.ticks.x = element_blank()
  )

ggsave(paste0(output_plot_prefix, ".qv.png"), plot = p_qv, width = plot_width, height = plot_height, limitsize = FALSE)

#barplots
#bars of regions ordered by tpr pct
#make summary
tpr_summary <- tpr_summary %>%
  mutate(accuracy = case_when(
    tpr_pct >= 95 ~ "high (>= 95%)",
    tpr_pct >= 80 ~ "mid (>= 80%)",
    TRUE ~ "low (< 80%)"
  ))

#barplot
#sort by tpr

tpr_summary_sorted <- tpr_summary %>%
  arrange(desc(tpr_pct))

#calculate number of rows and bars per row
max_bars_per_row <- 100
num_regions_tpr <- nrow(tpr_summary_sorted)
num_rows_tpr <- ceiling(num_regions_tpr / max_bars_per_row)
#how many bars per row
bars_per_row <- ceiling(num_regions_tpr / num_rows_tpr)

# Calculate the total number of samples for each region for TPR
tpr_region_samples <- data_long %>%
  select(region, sample.id) %>%
  distinct() %>%
  count(region, name = "sample_count")

tpr_bar_plots <- list()

#split plot into lines
for (i in 1:num_rows_tpr) {
  start_idx <- (i-1) * bars_per_row + 1
  end_idx <- min(i * bars_per_row, num_regions_tpr)
  if (start_idx > num_regions_tpr) break
  row_data <- tpr_summary_sorted[start_idx:end_idx, ]
  
  # Get the sample counts for these regions
  region_counts_tpr <- tpr_region_samples %>%
    filter(region %in% row_data$region)
  
  # Create a mapping dataframe that has both region names and tpr values
  region_mapping <- row_data %>%
    select(region, tpr_pct)
  
  # Join the mapping to get tpr_pct values in the region_counts_tpr dataframe
  region_counts_with_tpr <- region_counts_tpr %>%
    left_join(region_mapping, by = "region")
  
  p <- ggplot(row_data, aes(x = reorder(region, tpr_pct, decreasing = TRUE), y = tpr_pct, fill = accuracy)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("high (>= 95%)" = "#4CAF50", "mid (>= 80%)" = "#FFC107", "low (< 80%)" = "#F44336"),
      breaks = c("high (>= 95%)", "mid (>= 80%)", "low (< 80%)")
    ) +
    # Display sample counts vertically
    geom_text(
      data = region_counts_with_tpr,
      aes(x = region, y = 101, label = sample_count),
      inherit.aes = FALSE,
      angle = 90,  # Make text vertical
      hjust = -0.1,   # Align to bottom of text
      vjust = 0.5, # Center horizontally
      size = 5.8
    ) +
    labs(
      x = if(i == num_rows_tpr) "region" else "",
      y = "tpr (%)"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = 18),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.spacing.x = unit(2, "cm"),
      legend.key.size = unit(0.8, "line"),
      legend.text = element_text(size = 10,margin = margin(r = 10))
    ) +
    # Make room for the vertical labels
    scale_y_continuous(
      limits = c(0, 120),
      breaks = seq(0, 100, by = 10)
    )
  tpr_bar_plots[[i]] <- p
}

tpr_combined_plot <- plot_grid(plotlist = tpr_bar_plots, ncol = 1, align = 'v', axis = 'lr')

#plot dimensions
tpr_plot_width <- max(15, min(30, bars_per_row * 0.25))
tpr_plot_height <- 4.9 * num_rows_tpr

ggsave(paste0(output_plot_prefix, ".tpr_bar.png"), plot = tpr_combined_plot,width = tpr_plot_width, height = tpr_plot_height, limitsize = FALSE)

#qv_bar
#this should be similar to what locityper does in the paper
#the same logic above applies below

qv_data <- data_long %>%
  filter(metric %in% c("qv.1", "qv.2")) %>%
  select(region, sample.id, metric, metric.values)

qv_data <- qv_data %>%
  mutate(
    quality = factor(
      case_when(
        metric.values > 33            ~ "high: >33",
        metric.values > 23            ~ "mid: >23, <=33",
        metric.values > 17            ~ "low: >17, <= 23",
        metric.values <= 17           ~ "very low: <= 17",
        TRUE                          ~ "unknown"
      ),
      levels = c("unknown","very low: <= 17", "low: >17, <= 23","mid: >23, <=33", "high: >33")
    )
  )

qv_summary <- qv_data %>%
  count(region, quality) %>%
  group_by(region) %>%
  mutate(
    total = sum(n),
    percent = n / total * 100
  ) %>%
  ungroup()

# Calculate the total number of samples for each region
# Divide by 2 because we have qv.1 and qv.2 for each sample
region_samples <- qv_data %>%
  select(region, sample.id) %>%
  distinct() %>%
  count(region, name = "sample_count")

region_order_df <- qv_summary %>%
  filter(quality == "high: >33") %>%
  group_by(region) %>%
  summarize(high_pct = sum(percent)) %>%
  arrange(desc(high_pct))

all_regions <- unique(qv_summary$region)
missing_regions <- setdiff(all_regions, region_order_df$region)
if (length(missing_regions) > 0) {
  missing_df <- data.frame(region = missing_regions, high_pct = 0)
  region_order_df <- rbind(region_order_df, missing_df)
}

region_order <- region_order_df$region

qv_summary_sorted <- qv_summary %>%
  mutate(region = factor(region, levels = region_order)) %>%
  arrange(region)

num_regions_qv <- length(unique(qv_summary_sorted$region))
num_rows_qv <- ceiling(num_regions_qv / max_bars_per_row)
qv_bars_per_row <- ceiling(num_regions_qv / num_rows_qv)

qv_bar_plots <- list()
for (i in 1:num_rows_qv) {
  start_idx <- (i-1) * qv_bars_per_row + 1
  end_idx <- min(i * qv_bars_per_row, num_regions_qv)
  if (start_idx > num_regions_qv) break
  regions_in_row <- region_order[start_idx:end_idx]
  
  row_data <- qv_summary_sorted %>%
    filter(region %in% regions_in_row) %>%
    mutate(region = factor(region, levels = regions_in_row)) # Reorder within this subset
  
  # Get the sample counts for these regions
  region_counts <- region_samples %>%
    filter(region %in% regions_in_row)
  
  p <- ggplot(row_data, aes(x = region, y = percent, fill = quality)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c(
        "high: >33" = "#4CAF50",
        "mid: >23, <=33" = "#FFC107",
        "low: >17, <= 23" = "#FF8C00",
        "very low: <= 17" = "#F44336"
      )
    ) +
    # Display sample counts vertically
    geom_text(
      data = region_counts,
      aes(x = region, y = 101, label = sample_count),
      inherit.aes = FALSE,
      angle = 90,  # Make text vertical
      hjust = -0.1,   # Align to bottom of text
      vjust = 0.5, # Center horizontally
      size = 5.8
    ) +
    labs(
      x = if(i == num_rows_qv) "region" else "", 
      y = "pct of QVs"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text = element_text(size = 18),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.spacing.x = unit(2, "cm"),
      legend.text = element_text(size = 10,margin = margin(r = 10)),
      legend.key.size = unit(0.8, "line")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    # Make room for the vertical labels
    scale_y_continuous(
      limits = c(0, 120),
      breaks = seq(0, 100, by = 10)
    )
  qv_bar_plots[[i]] <- p
}

#calculate dimensions
qv_plot_width <- max(15, min(30, qv_bars_per_row * 0.25))
qv_plot_height <- 4.9 * num_rows_qv

#combine
qv_combined_plot <- plot_grid(plotlist = qv_bar_plots, ncol = 1, align = 'v', axis = 'lr')

ggsave(paste0(output_plot_prefix, ".qv_bar.png"), plot = qv_combined_plot, width = qv_plot_width, height = qv_plot_height, limitsize = FALSE)
#done