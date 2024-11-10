#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(ggplot2)

# Set data.table threads
setDTthreads(1)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_time_plot <- args[2]
output_mem_plot <- args[3]

# Helper functions
get_method <- function(x) {
  tail(unlist(strsplit(x, ".", fixed = TRUE)), 3)[1]
}

get_region <- function(x) {
  y <- tail(unlist(strsplit(x, ".", fixed = TRUE)), 4)[1]
  tail(unlist(strsplit(y, "/")), 1)
}

# List and categorize input files
files <- sort(list.files(input_dir, pattern = ".benchmark.txt", full.names = TRUE))
files_wregion <- grep("^[^_]+_[^_]+_[^_]+_[^_]+[^.]_", files, value = TRUE)
files_noregion <- grep("^[^_]+_[^_]+_[^_]+_[^_]+[^.]_", files, invert = TRUE, value = TRUE)

# Process files with regions
method <- unique(sapply(files_wregion, get_method))
region <- unique(sapply(files_wregion, get_region))
benchmark_list <- list()

for (m in method) {
  for (r in region) {
    f1 <- files_wregion[grep(m, files_wregion)]
    f2 <- f1[grep(r, f1)]
    if (length(f2) != 0) {
      df <- rbindlist(lapply(f2, fread))
      benchmark_list[[length(benchmark_list) + 1]] <- data.frame(
        region = r,
        method = m,
        value.time = df$s,
        value.mem = df$max_rss
      )
    }
  }
}

# Process files without regions
for (f in files_noregion) {
  method <- tail(unlist(strsplit(basename(f), "\\.")), 4)[1]
  df <- fread(f)
  benchmark_list[[length(benchmark_list) + 1]] <- data.frame(
    region = "Shared",
    method = method,
    value.time = df$s,
    value.mem = df$max_rss
  )
}

benchmark_df <- rbindlist(benchmark_list)

# Create and save runtime plot
p_time <- ggplot(benchmark_df, aes(x = method, y = value.time, fill = region)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000)) +
  labs(x = "Step", y = "Runtime (s)")

ggsave(output_time_plot, plot = p_time, width = 20)

# Create and save memory usage plot
p_mem <- ggplot(benchmark_df, aes(x = method, y = value.mem, fill = region)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "top") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  labs(x = "Step", y = "Max RSS (MB)")

ggsave(output_mem_plot, plot = p_mem, width = 20)