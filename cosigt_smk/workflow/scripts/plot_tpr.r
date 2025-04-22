#!/usr/bin/Rscript

library(ggplot2)
library(ggrepel)
library(data.table)
library(scales)
library(dplyr)
library(tidyr)
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
                             tpr_summary$tp_count, 
                             tpr_summary$total_count, 
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
# Calculate number of columns (max 20 per row)
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

#bars of regions ordered by tpr pct
#make summary
tpr_summary <- tpr_summary %>%
  mutate(performance_category = case_when(
    tpr_pct >= 95 ~ "high (>= 95%)",
    tpr_pct >= 80 ~ "mid (>= 80%)",
    TRUE ~ "low (< 80%)"
  ))

#barplot
p_bar_tpr <- ggplot(tpr_summary, aes(x = reorder(region, tpr_pct, decreasing = T), y = tpr_pct, fill = performance_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("high (>= 95%)" = "#4CAF50", "mid (>= 80%)" = "#FFC107", "low (< 80%)" = "#F44336"),
    breaks = c("high (>= 95%)", "mid (>= 80%)", "low (< 80%)")
  ) +
  labs(
    x = "region",
    y = "tpr (%)",
    fill = "performance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 5)
  )

ggsave(paste0(output_plot_prefix, ".tpr_bar.png"), plot = p_bar_tpr, width = plot_width, height = plot_height, limitsize = FALSE)

#qv_bar
#this should be similar to what locityper does in the paper
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

region_order <- qv_summary %>%
  filter(quality == "high: >33") %>%
  arrange(desc(percent)) %>%
  pull(region)

qv_summary$region <- factor(qv_summary$region, levels = region_order)

p_qv_bar<-ggplot(qv_summary, aes(x = region, y = percent, fill = quality)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "high: >33" = "#4CAF50",
      "mid: >23, <=33" = "#FFC107",
      "low: >17, <= 23" = "#FF8C00",
      "very low: <= 17" = "#F44336"
    )
  ) +
  labs(
    x = "region",
    y = "pct of QVs",
    fill = "accuracy",
    title = "qv metric - quality distribution"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "top"
  )

ggsave(paste0(output_plot_prefix, ".qv_bar.png"), plot = p_qv_bar, width = plot_width, height = plot_height, limitsize = FALSE)
#done