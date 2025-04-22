#!/usr/bin/env Rscript

library(data.table)
library(dbscan)
library(ggplot2)
library(ggrepel)

# Set data.table threads
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
bed_in <- args[1]
bed_out <- args[2]
plot_out <- sub("\\.bed$", "_length_distribution.png", bed_out)  # Create plot filename based on output

df <- fread(bed_in)
df$length <- df$V3 - df$V2

# create distance matrix directly from lengths
create_length_distance_matrix <- function(lengths) {
  n <- length(lengths)
  dist_matrix <- matrix(0, nrow=n, ncol=n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      pct_diff <- 1 - min(lengths[i], lengths[j]) / max(lengths[i], lengths[j])
      dist_matrix[i, j] <- pct_diff
      dist_matrix[j, i] <- pct_diff
    }
  }
  
  return(as.dist(dist_matrix))
}

dist_matrix <- create_length_distance_matrix(df$length)

#log-scale proportion (10% usually) is n haps increases
calculate_min_pts <- function(n, base_proportion=0.1, min_value=1) {
  if (n <= 20) {
    return(min_value)
  } else if (n <= 200) {
    return(ceiling(n * base_proportion))
  } else {
    return(ceiling(log10(n) * 3))
  }
}

min_pts <- calculate_min_pts(nrow(df))
res <- dbscan(dist_matrix, eps=0.2, minPts=min_pts)$cluster
outlier_indices <- which(res == 0)

plot_df <- copy(df)
plot_df$outlier <- FALSE
if (length(outlier_indices) > 0) {
  plot_df$outlier[outlier_indices] <- TRUE
}

set.seed(42)
plot_df$x_jitter <- jitter(rep(1, nrow(plot_df)), amount=0.1)

#viz
p <- ggplot(plot_df, aes(x=x_jitter, y=length)) +
  geom_violin(aes(x=x_jitter,group=1), fill="lightblue", alpha=0.5) +
  geom_boxplot(aes(x=x_jitter,group=1), width=0.2, alpha=0.7) +
  geom_point(aes(color=outlier), alpha=0.7, size=3) +
  scale_color_manual(values=c("FALSE"="darkblue", "TRUE"="red")) +
  labs(title="distibution of region lengths",
       subtitle=paste0("outliers detected: ", sum(plot_df$outlier), "/", nrow(plot_df), 
                      " (minPts=", min_pts, ", eps=0.2)"),
       y="length (bp)", 
       x="") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

if (sum(plot_df$outlier) > 0) {
  outlier_df <- plot_df[plot_df$outlier == TRUE,]
  p <- p + geom_text_repel(
    data=outlier_df,
    aes(label=V1),
    size=3,
    box.padding = 0.5,
    point.padding = 0.5,
    force=2,
    segment.color="red",
    min.segment.length = 0.1
  )
}

ggsave(plot_out, p, width=8, height=10)
if (length(outlier_indices) > 0) {
  df <- df[-outlier_indices]
}
fwrite(df, bed_out, row.names = FALSE, col.names = FALSE, sep = "\t")
