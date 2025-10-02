#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

table_in <- args[1]
output_prefix <- args[2]
node_length_file <- args[3]
panplexity_mask_file <- args[4]

library(data.table)
library(ggplot2)
library(gridExtra)

df <- fread(table_in, sep = "\t", header = TRUE, check.names = FALSE)
node_lengths <- fread(node_length_file, sep = "\t", header = FALSE, col.names = c("node", "length"))

# Panplexity mask
panplexity_mask <- fread(panplexity_mask_file, sep="\t", header=FALSE, col.names=c("panplexity_mask"))
panplexity_mask$node<-node_lengths$node

# Coverage stats by node (mean)
df_long <- melt(df, id.vars = "path.name", variable.name = "node", value.name = "coverage")
coverage_stats <- df_long[, .(mean_coverage = mean(coverage, na.rm = TRUE)), by = node]

# Merge length and mask
node_stats <- merge(coverage_stats, node_lengths, by="node")
node_stats <- merge(node_stats, panplexity_mask, by="node")
node_stats[, node := as.numeric(sub("node\\.", "", node))]
setorder(node_stats, node)
node_stats[, node_factor := factor(node, levels = unique(node[order(node)]))]

# IQR-based coverage mask (1 = include, 0 = mask)
Q1 <- quantile(node_stats$mean_coverage, 0.25, na.rm = TRUE)
Q3 <- quantile(node_stats$mean_coverage, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val
node_stats[, coverage_mask := as.integer(mean_coverage >= lower_bound & mean_coverage <= upper_bound)]

# Final mask: only 1 if both masks are 1 (LOGICAL AND)
node_stats[, final_mask := as.integer(panplexity_mask == 1 & coverage_mask == 1)]

node_stats[, mask_type :=
  fifelse(panplexity_mask == 1 & coverage_mask == 1, "unmasked",
    fifelse(panplexity_mask == 1 & coverage_mask == 0, "coverage_masked",
      fifelse(panplexity_mask == 0 & coverage_mask == 1, "panplexity_masked", "both_masked")
    )
  )
]
node_stats$mask_type <- factor(node_stats$mask_type, levels = c("unmasked", "panplexity_masked", "coverage_masked", "both_masked"))

# Plots
p1 <- ggplot(node_stats, aes(x = length, fill = mask_type)) +
  geom_density(alpha = 0.6) +
  labs(title = "Node Length Density", x = "Length", y = "Density", fill = "Mask type") +
  theme_minimal()

p2 <- ggplot(node_stats, aes(x = mean_coverage, fill = mask_type)) +
  geom_density(alpha = 0.6) +
  labs(title = "Mean Coverage Density", x = "Mean coverage", y = "Density", fill = "Mask type") +
  theme_minimal()

p3 <- ggplot(node_stats, aes(x = node, y = mask_type, color = mask_type)) +
  geom_point(size = 1.25, alpha=0.7) +
  labs(title = "Masking Status by Node Index", x = "Node index", y = "Mask type", color = "Mask type") +
  theme_minimal()

comb_plot <- grid.arrange(p1, p2, p3, ncol = 1)
ggsave(paste0(output_prefix, "_mask_distributions.png"), comb_plot, width = 10, height = 20, dpi = 400)

mask_table <- node_stats[, .(node, panplexity_mask, coverage_mask, final_mask)]
fwrite(mask_table, paste0(output_prefix, "_node_masks.tsv"), sep = "\t", quote = FALSE)

total_length <- sum(node_stats$length, na.rm = TRUE)
summarize_mask <- function(mask_col) sum(node_stats$length[mask_col == 0], na.rm = TRUE)
panp_masked <- summarize_mask(node_stats$panplexity_mask)
cov_masked <- summarize_mask(node_stats$coverage_mask)
both_masked <- summarize_mask(node_stats$final_mask)
unmasked <- sum(node_stats$length[node_stats$panplexity_mask == 1 & node_stats$coverage_mask == 1], na.rm = TRUE)

summary_tab <- data.table(
  mask = c("panplexity_masked", "coverage_masked", "both_masked", "unmasked"),
  bases = c(panp_masked, cov_masked, both_masked, unmasked),
  pct_bases = round(100 * c(panp_masked, cov_masked, both_masked, unmasked) / total_length, 2)
)
fwrite(summary_tab, paste0(output_prefix, "_mask_summary.tsv"), sep = "\t", quote = FALSE)

fwrite(as.list(node_stats$final_mask),
       paste0(output_prefix, ".mask.tsv"),
       sep = "\n", col.names = FALSE)