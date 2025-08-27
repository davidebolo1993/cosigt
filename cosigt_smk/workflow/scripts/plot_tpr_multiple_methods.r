#!/usr/bin/Rscript

# USAGE
# Concatenate "tpr_qv.tsv" files for each set
# cat $dir_base/linguistic-k16-t1/benchmark/*/*/tpr_qv.tsv > linguistic-k16-t1.tsv
# cat $dir_base/linguistic-k16-tauto/benchmark/*/*/tpr_qv.tsv > linguistic-k16-tauto.tsv
# Then run
# Rscript plot_tpr_multiple_methods.R prefix-you-like regions linguistic-k16-t1.tsv linguistic-k16-tauto.tsv

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
output_plot_prefix <- args[1]
annot_bed <- args[2]
annot_tsv <- fread(annot_bed)
annot_tsv <- annot_tsv %>%
  mutate(id = paste0(as.character(V1), "_", as.character(V2), "_", as.character(V3)))

# All remaining arguments are TSV files
file_paths <- args[3:length(args)]
n_lists <- length(file_paths)

if (n_lists < 1) {
  stop("At least one TSV file must be provided")
}

# Extract method names from file paths
list_names <- character(n_lists)
for (i in 1:n_lists) {
  # Get basename without extension
  filename <- basename(file_paths[i])
  # Remove .tsv extension
  list_names[i] <- sub("\\.tsv$", "", filename)
  # Clean up common patterns (optional - customize as needed)
  list_names[i] <- gsub("_", "-", list_names[i])
  list_names[i] <- gsub("tpr-qv", "", list_names[i])
  list_names[i] <- gsub("^-|-$", "", list_names[i])  # Remove leading/trailing dashes
  # If empty after cleaning, use a default name
  if (list_names[i] == "" || is.na(list_names[i])) {
    list_names[i] <- paste0("Method", i)
  }
}

print(paste("Detected methods:", paste(list_names, collapse=", ")))

# Read TPR files for each method
all_tpr_dfs <- list()
for (i in 1:n_lists) {
  print(paste("Reading", file_paths[i], "..."))
  
  # Read the file, handling potential duplicate headers from concatenation
  # First read to detect column names
  first_line <- readLines(file_paths[i], n = 1)
  
  # Read the entire file
  tpr_df <- fread(file_paths[i], header = TRUE)
  
  # Remove any rows that are duplicate headers (they'll have column names as values)
  # Assuming first column is 'sample.id' or similar
  col_names <- names(tpr_df)
  if ("sample.id" %in% col_names) {
    tpr_df <- tpr_df[sample.id != "sample.id"]
  } else if (length(col_names) > 0) {
    # Generic approach: remove rows where first column equals the column name
    first_col <- col_names[1]
    tpr_df <- tpr_df[get(first_col) != first_col]
  }
  
  # Convert columns to appropriate types (they might be character after removing headers)
  numeric_cols <- c("edr.1", "edr.2", "qv.1", "qv.2", "n.clust")
  for (col in numeric_cols) {
    if (col %in% names(tpr_df)) {
      tpr_df[[col]] <- as.numeric(tpr_df[[col]])
    }
  }
  
  # Add method identifier
  tpr_df$list_id <- list_names[i]
  
  # Remove any rows with NA in critical columns
  tpr_df <- tpr_df[!is.na(sample.id) & !is.na(region)]
  
  all_tpr_dfs[[i]] <- tpr_df
  print(paste("Loaded", nrow(tpr_df), "valid rows from", list_names[i]))
}

# Combine all methods
combined_tpr_df <- do.call(rbind, all_tpr_dfs)

# Convert to long format
data_long <- combined_tpr_df %>%
  select(sample.id, region, n.clust, edr.1, edr.2, qv.1, qv.2, tpr, list_id) %>% 
  pivot_longer(cols = c(edr.1, edr.2, qv.1, qv.2), names_to = "metric", values_to = "metric.values")

# Add annotations if available
data_long$pos <- data_long$region
if (any(annot_tsv$V4 != "unknown")) {
  # Match and update regions
  matched_regions <- annot_tsv$V4[match(data_long$pos, annot_tsv$id)]
  # Only update if match found (keep original if no match)
  data_long$region <- ifelse(is.na(matched_regions), data_long$region, matched_regions)
}

# Remove any rows where region is NA or empty
data_long <- data_long %>%
  filter(!is.na(region) & region != "" & region != "NA")

print(paste("Unique regions found:", paste(unique(data_long$region), collapse=", ")))

# Filter for QV data only
qv_data <- data_long %>%
  filter(metric %in% c("qv.1", "qv.2")) %>%
  select(region, sample.id, metric, metric.values, list_id) %>%
  filter(!is.na(metric.values))  # Remove NA values

# Categorize QV values
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

# Create summary statistics per method
qv_summary <- qv_data %>%
  count(region, quality, list_id) %>%
  group_by(region, list_id) %>%
  mutate(
    total = sum(n),
    percent = n / total * 100
  ) %>%
  ungroup()

# Count samples per region per method
region_samples <- qv_data %>%
  select(region, sample.id, list_id) %>%
  distinct() %>%
  count(region, list_id, name = "sample_count")

# Order regions by average percentage of high quality QVs across all methods
region_order_df <- qv_summary %>%
  filter(quality == "high: >33") %>%
  group_by(region) %>%
  summarize(avg_high_pct = mean(percent, na.rm = TRUE)) %>%
  arrange(desc(avg_high_pct))

# Add regions with no high quality QVs
all_regions <- unique(qv_summary$region)
missing_regions <- setdiff(all_regions, region_order_df$region)
if (length(missing_regions) > 0) {
  missing_df <- data.frame(region = missing_regions, avg_high_pct = 0)
  region_order_df <- rbind(region_order_df, missing_df)
}

region_order <- region_order_df$region

# Sort summary by region order
qv_summary_sorted <- qv_summary %>%
  mutate(region = factor(region, levels = region_order)) %>%
  arrange(region, list_id)

# Save summary table
fwrite(qv_summary_sorted, file=paste0(output_plot_prefix, ".qv_bar.tsv"), 
       sep="\t", quote=F, col.names=T, row.names=F)

# Plot parameters
max_regions_per_row <- 20  # Reduced since we have multiple bars per region
num_regions <- length(unique(qv_summary_sorted$region))
num_rows <- ceiling(num_regions / max_regions_per_row)
regions_per_row <- ceiling(num_regions / num_rows)

# Create plots for each row
qv_bar_plots <- list()
for (i in 1:num_rows) {
  start_idx <- (i-1) * regions_per_row + 1
  end_idx <- min(i * regions_per_row, num_regions)
  if (start_idx > num_regions) break
  regions_in_row <- region_order[start_idx:end_idx]
  
  row_data <- qv_summary_sorted %>%
    filter(region %in% regions_in_row) %>%
    mutate(region = factor(region, levels = regions_in_row))
  
  region_counts <- region_samples %>%
    filter(region %in% regions_in_row)
  
  # Create grouped bar plot
  p <- ggplot(row_data, aes(x = list_id, y = percent, fill = quality)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ region, ncol = length(regions_in_row), scales = "free_x") +
    scale_fill_manual(
      values = c(
        "high: >33" = "#4CAF50",
        "mid: >23, <=33" = "#FFC107",
        "low: >17, <= 23" = "#FF8C00",
        "very low: <= 17" = "#F44336"
      )
    ) +
    # Add sample counts
    geom_text(
      data = region_counts,
      aes(x = list_id, y = 101, label = sample_count),
      inherit.aes = FALSE,
      angle = 45,
      hjust = 0,
      vjust = 0.5,
      size = 3.5
    ) +
    labs(
      x = if(i == num_rows) "Method" else "", 
      y = "pct of QVs",
      fill = "Quality"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.key.size = unit(0.8, "line"),
      panel.spacing = unit(0.5, "lines")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_y_continuous(
      limits = c(0, 120),
      breaks = seq(0, 100, by = 25)
    )
  
  qv_bar_plots[[i]] <- p
}

# Calculate plot dimensions
plot_width <- max(12, min(30, regions_per_row * n_lists * 0.8))
plot_height <- 6 * num_rows

# Combine and save plot
if (length(qv_bar_plots) == 1) {
  qv_combined_plot <- qv_bar_plots[[1]]
} else {
  qv_combined_plot <- plot_grid(plotlist = qv_bar_plots, ncol = 1, align = 'v', axis = 'lr')
}

ggsave(paste0(output_plot_prefix, ".qv_bar.png"), plot = qv_combined_plot, 
       width = plot_width, height = plot_height, limitsize = FALSE)

# Create comparison plot showing only high QV percentages
high_qv_comparison <- qv_summary %>%
  filter(quality == "high: >33") %>%
  select(region, list_id, percent) %>%
  complete(region, list_id, fill = list(percent = 0))

# Order methods consistently
list_order <- unique(list_names)
high_qv_comparison$list_id <- factor(high_qv_comparison$list_id, levels = list_order)

comparison_plot <- ggplot(high_qv_comparison, 
                          aes(x = reorder(region, -percent), y = percent, fill = list_id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Region",
    y = "Percentage of High QV (>33)",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  ) +
  scale_fill_brewer(palette = "Set1")

ggsave(paste0(output_plot_prefix, ".qv_comparison.png"), 
       plot = comparison_plot, 
       width = max(12, min(30, length(unique(high_qv_comparison$region)) * 0.3)), 
       height = 8)

# NEW: Simple summary plot - count of regions with >=80% high QV
# Calculate percentage of high QV for each region and method
high_qv_by_region <- qv_summary %>%
  filter(quality == "high: >33") %>%
  select(region, list_id, percent) %>%
  complete(region, list_id, fill = list(percent = 0))

# Count regions meeting the 80% threshold for each method
regions_above_threshold <- high_qv_by_region %>%
  filter(percent >= 80) %>%
  count(list_id, name = "n_regions_high_qv") %>%
  complete(list_id = unique(high_qv_by_region$list_id), fill = list(n_regions_high_qv = 0))

# Add percentage for context
total_regions <- length(unique(high_qv_by_region$region))
regions_above_threshold <- regions_above_threshold %>%
  mutate(
    percent_regions = (n_regions_high_qv / total_regions) * 100,
    label = paste0(n_regions_high_qv, "\n(", round(percent_regions, 1), "%)")
  )

# Order methods consistently
regions_above_threshold$list_id <- factor(regions_above_threshold$list_id, levels = list_order)

# Create the simple bar plot
summary_plot <- ggplot(regions_above_threshold, 
                       aes(x = list_id, y = n_regions_high_qv, fill = list_id)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), vjust = -0.5, size = 4) +
  labs(
    x = "Method",
    y = "Number of regions with ≥80% high QV (>33)",
    title = paste0("Regions with high-quality assemblies (≥80% QV>33)\nTotal regions: ", total_regions)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none"  # Remove legend since colors are redundant with x-axis
  ) +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, max(regions_above_threshold$n_regions_high_qv) * 1.15)  # Add space for labels

ggsave(paste0(output_plot_prefix, ".qv_summary.png"), 
       plot = summary_plot, 
       width = max(8, n_lists * 2), 
       height = 8)

# Save the summary statistics to a file
fwrite(regions_above_threshold, 
       file=paste0(output_plot_prefix, ".qv_summary_stats.tsv"), 
       sep="\t", quote=F, col.names=T, row.names=F)

print("Summary:")
print(paste("Total regions:", num_regions))
for (i in 1:n_lists) {
  n_samples <- length(unique(qv_data[qv_data$list_id == list_names[i], ]$sample.id))
  n_high_qv <- regions_above_threshold$n_regions_high_qv[regions_above_threshold$list_id == list_names[i]]
  pct_high_qv <- regions_above_threshold$percent_regions[regions_above_threshold$list_id == list_names[i]]
  print(paste(" -", list_names[i], ":", n_samples, "samples,", 
              n_high_qv, "regions with ≥80% high QV",
              paste0("(", round(pct_high_qv, 1), "%)")))
}
print("Plots saved:")
print(paste(" -", paste0(output_plot_prefix, ".qv_bar.png"), "(detailed QV distribution)"))
print(paste(" -", paste0(output_plot_prefix, ".qv_comparison.png"), "(high QV comparison)"))
print(paste(" -", paste0(output_plot_prefix, ".qv_summary.png"), "(regions meeting 80% threshold)"))
print(paste(" -", paste0(output_plot_prefix, ".qv_summary_stats.tsv"), "(summary statistics)"))