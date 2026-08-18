#!/usr/bin/env Rscript
# Per-gene summaries of the benchmark QV table, following the cosigt paper's
# supplementary figures: a stacked bar per gene showing how the per-haplotype
# values distribute across bands, genes ordered by the share in the best band,
# and the number of samples printed above each bar.
#
#   Rscript plot_qv.r <benchmark.qv.tsv> <output.png> <mode> [max_bars_per_row]
#
# One plot, chosen by the benchmark mode that produced the table:
#   leave_zero_out  QV_pred banded into four quality bands
#   leave_all_out   QVfrac binned into quintiles, because
#                   there the achieved QV alone cannot separate genotyping error
#                   from a panel that simply held nothing closer

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_qv.r <benchmark.qv.tsv> <output.png> <mode> [max_bars_per_row]")
}

qv_table <- args[[1]]
output_png <- args[[2]]
mode <- args[[3]]
max_bars_per_row <- if (length(args) >= 4) as.numeric(args[[4]]) else 30

if (!mode %in% c("leave_zero_out", "leave_all_out")) {
  stop("mode must be leave_zero_out or leave_all_out, got: ", mode)
}

# QV bands, colour-blind safe (Tol)
qv_levels <- c("very low: <= 17", "low: >17, <= 23", "mid: >23, <=33", "high: >33")
qv_colors <- c(
  "very low: <= 17" = "#DDCC77",
  "low: >17, <= 23" = "#999933",
  "mid: >23, <=33"  = "#117733",
  "high: >33"       = "#44AA99"
)

# QVfrac quintiles, viridis mako. Hardcoded so the plot does not depend on the
# viridis package being present in the container; produced by
# rev(viridis(n = 5, option = "mako", direction = -1, begin = 0.2, end = 0.9)).
frac_levels <- c("Q1: 0.0-0.2", "Q2: 0.2-0.4", "Q3: 0.4-0.6", "Q4: 0.6-0.8", "Q5: 0.8-1.0")
frac_colors <- c(
  "Q5: 0.8-1.0" = "#382A54",
  "Q4: 0.6-0.8" = "#3B5698",
  "Q3: 0.4-0.6" = "#3488A6",
  "Q2: 0.2-0.4" = "#43BAAD",
  "Q1: 0.0-0.2" = "#A8E1BC"
)

placeholder <- function(path, message) {
  ggsave(
    path,
    ggplot() + annotate("text", x = 0, y = 0, label = message, size = 6) + theme_void(),
    width = 10, height = 5, limitsize = FALSE
  )
}

# Shared layout: stacked percentage bars per gene, ordered by the share in
# `best_level`, sample counts above each bar, split into rows of at most
# max_bars_per_row via facets so the legend stays shared.
banded_bar_plot <- function(dt, band_col, levels_vec, colors_vec, best_level,
                            legend_title, y_title, path) {
  dt <- dt[!is.na(get(band_col))]
  if (nrow(dt) == 0) {
    placeholder(path, "No values available")
    return(invisible(NULL))
  }

  summary_dt <- dt[, .N, by = c("label", band_col)]
  summary_dt[, percent := 100 * N / sum(N), by = label]

  best <- summary_dt[get(band_col) == best_level, .(label, pct_best = percent)]
  missing <- setdiff(unique(summary_dt$label), best$label)
  if (length(missing) > 0) {
    best <- rbind(best, data.table(label = missing, pct_best = 0))
  }
  gene_order <- best[order(-pct_best)]$label
  summary_dt[, label := factor(label, levels = gene_order)]

  counts <- dt[, .(n_samples = uniqueN(sample)), by = label]
  counts[, label := factor(label, levels = gene_order)]

  n_genes <- length(gene_order)
  n_rows <- max(1, ceiling(n_genes / max_bars_per_row))
  bars_per_row <- ceiling(n_genes / n_rows)
  row_of <- data.table(
    label = factor(gene_order, levels = gene_order),
    row_id = ((seq_len(n_genes) - 1) %/% bars_per_row) + 1
  )
  summary_dt <- merge(summary_dt, row_of, by = "label")
  counts <- merge(counts, row_of, by = "label")

  p <- ggplot(summary_dt, aes(x = label, y = percent, fill = .data[[band_col]])) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    geom_text(
      data = counts,
      aes(x = label, y = 101, label = n_samples),
      inherit.aes = FALSE,
      angle = 90, hjust = -0.1, vjust = 0.5, size = 3
    ) +
    scale_fill_manual(values = colors_vec, name = legend_title, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 115),
      breaks = seq(0, 100, by = 25),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Gene", y = y_title) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 1))

  if (n_rows > 1) {
    p <- p + facet_wrap(~ row_id, ncol = 1, scales = "free_x")
  }

  ggsave(
    path, p,
    width = max(8, bars_per_row * 0.28 + 3),
    height = max(4.5, n_rows * 4),
    dpi = 300, limitsize = FALSE
  )
}

df <- fread(qv_table, header = TRUE)
df <- df[!is.na(QV_1_pred) & !is.na(QV_2_pred)]

if (nrow(df) == 0) {
  placeholder(output_png, "No samples with true haplotypes available for comparison")
  quit(status = 0)
}

# gene_name is column 4 of the regions BED; fall back to the region id so that
# unannotated regions stay distinct rather than collapsing into one bar
df[, label := ifelse(
  is.na(gene_name) | gene_name %in% c("", "unknown", "NA"),
  region,
  gene_name
)]

if (mode == "leave_zero_out") {
  # QV_pred banded into quality bands, one row per predicted haplotype
  long <- melt(
    df,
    id.vars = c("sample", "region", "label"),
    measure.vars = c("QV_1_pred", "QV_2_pred"),
    variable.name = "which",
    value.name = "QV"
  )
  long <- long[!is.na(QV)]
  long[, quality := factor(
    fifelse(QV > 33, "high: >33",
    fifelse(QV > 23, "mid: >23, <=33",
    fifelse(QV > 17, "low: >17, <= 23", "very low: <= 17"))),
    levels = qv_levels
  )]

  banded_bar_plot(
    long, "quality", qv_levels, qv_colors, "high: >33",
    expression(QV[pred]), expression("% of" ~ QV[pred]), output_png
  )
} else {
  # QVfrac in quintiles, one row per sample
  if (!"QVfrac" %in% names(df) || all(is.na(df$QVfrac))) {
    placeholder(output_png, "No QVfrac values: leave_all_out produced no comparable samples")
    quit(status = 0)
  }
  frac <- df[!is.na(QVfrac)]
  frac[, quintile := cut(
    round(QVfrac, 2),
    breaks = seq(0, 1, 0.2),
    include.lowest = TRUE,
    labels = frac_levels
  )]

  banded_bar_plot(
    frac, "quintile", frac_levels, frac_colors, "Q5: 0.8-1.0",
    expression(QV[frac]), expression("% of" ~ QV[frac]), output_png
  )
}
