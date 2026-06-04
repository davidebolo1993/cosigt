#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(reshape2)
  library(rjson)
  library(dbscan)
  library(cluster)
})

setDTthreads(1)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript workflow/scripts/cluster.r input.tsv output.json [options]",
      "",
      "Legacy positional form is still accepted:",
      "  Rscript workflow/scripts/cluster.r input.tsv output.json similarity_threshold region_similarity levels [options]",
      "",
      "Options:",
      "  --similarity-threshold automatic|N Fixed similarity or automatic eps [automatic]",
      "  --levels N                      Recursive clustering depth [1]",
      "  --eps-selection quality|legacy  Automatic eps selection mode [quality]",
      "  --eps-min N                     Minimum eps considered for automatic selection [0]",
      "  --eps-max N                     Maximum eps considered for automatic selection [0.30]",
      "  --eps-step N                    Grid step for automatic eps scan [0.01]",
      "  --score-use-dbcv true|false     Use DBCV in automatic eps scoring [true]",
      "  --dbcv-dim N                    Dimensionality passed to DBCV for dist input [2]",
      "  --low-diversity-mpd-norm N      Allow one giant/one cluster below this mpd_norm [0.05]",
      "  --giant-cluster-fraction N      Penalize giant clusters above this fraction [0.85]",
      "  --small-cluster-size N          Size cutoff for small-cluster diagnostics [2]",
      "  --ambiguous-eps-ratio N         Small clusters this close to eps are ambiguous [1.25]",
      "  --plot-clusters true|false      Write a PCoA cluster plot [false]",
      "  --plot-label-max-cluster-size N Label haplotypes in clusters with <= N members [1]",
      "  --plot-width N                  Cluster plot width in inches [8]",
      "  --plot-height N                 Cluster plot height in inches [6]",
      "",
      sep = "\n"
    ),
    file = stderr()
  )
}

parse_bool <- function(value) {
  tolower(value) %in% c("true", "t", "1", "yes", "y")
}

parse_options <- function(extra_args) {
  opt <- list(
    eps_selection = "quality",
    similarity_threshold = "automatic",
    levels = 1,
    eps_min = 0,
    eps_max = 0.30,
    eps_step = 0.01,
    score_use_dbcv = TRUE,
    dbcv_dim = 2,
    low_diversity_mpd_norm = 0.05,
    giant_cluster_fraction = 0.85,
    small_cluster_size = 2,
    ambiguous_eps_ratio = 1.25,
    plot_clusters = FALSE,
    plot_label_max_cluster_size = 1,
    plot_width = 8,
    plot_height = 6
  )

  if (length(extra_args) == 0) {
    return(opt)
  }
  if (length(extra_args) %% 2 != 0) {
    usage()
    stop("Optional arguments must be supplied as --key value pairs.")
  }

  i <- 1
  while (i <= length(extra_args)) {
    key <- extra_args[[i]]
    value <- extra_args[[i + 1]]
    if (key == "--similarity-threshold") {
      opt$similarity_threshold <- value
    } else if (key == "--levels") {
      opt$levels <- as.integer(value)
    } else if (key == "--eps-selection") {
      opt$eps_selection <- value
    } else if (key == "--eps-min") {
      opt$eps_min <- as.numeric(value)
    } else if (key == "--eps-max") {
      opt$eps_max <- as.numeric(value)
    } else if (key == "--eps-step") {
      opt$eps_step <- as.numeric(value)
    } else if (key == "--score-use-dbcv") {
      opt$score_use_dbcv <- parse_bool(value)
    } else if (key == "--dbcv-dim") {
      opt$dbcv_dim <- as.numeric(value)
    } else if (key == "--low-diversity-mpd-norm") {
      opt$low_diversity_mpd_norm <- as.numeric(value)
    } else if (key == "--giant-cluster-fraction") {
      opt$giant_cluster_fraction <- as.numeric(value)
    } else if (key == "--small-cluster-size") {
      opt$small_cluster_size <- as.integer(value)
    } else if (key == "--ambiguous-eps-ratio") {
      opt$ambiguous_eps_ratio <- as.numeric(value)
    } else if (key == "--plot-clusters") {
      opt$plot_clusters <- parse_bool(value)
    } else if (key == "--plot-label-max-cluster-size") {
      opt$plot_label_max_cluster_size <- as.integer(value)
    } else if (key == "--plot-width") {
      opt$plot_width <- as.numeric(value)
    } else if (key == "--plot-height") {
      opt$plot_height <- as.numeric(value)
    } else {
      stop("Unknown option: ", key)
    }
    i <- i + 2
  }

  if (!opt$eps_selection %in% c("quality", "legacy")) {
    stop("--eps-selection must be 'quality' or 'legacy'.")
  }
  if (is.na(opt$levels) || opt$levels < 1) {
    stop("--levels must be a positive integer.")
  }
  if (!is.finite(opt$eps_min) || !is.finite(opt$eps_max) || opt$eps_min < 0 ||
      opt$eps_max < opt$eps_min) {
    stop("--eps-min/--eps-max must define a finite non-negative range.")
  }
  if (!is.finite(opt$eps_step) || opt$eps_step <= 0) {
    stop("--eps-step must be > 0.")
  }
  if (!is.finite(opt$dbcv_dim) || opt$dbcv_dim <= 0) {
    stop("--dbcv-dim must be > 0.")
  }
  if (!is.finite(opt$low_diversity_mpd_norm) || opt$low_diversity_mpd_norm < 0) {
    stop("--low-diversity-mpd-norm must be >= 0.")
  }
  if (!is.finite(opt$giant_cluster_fraction) || opt$giant_cluster_fraction <= 0 ||
      opt$giant_cluster_fraction >= 1) {
    stop("--giant-cluster-fraction must be between 0 and 1.")
  }
  if (is.na(opt$small_cluster_size) || opt$small_cluster_size < 1) {
    stop("--small-cluster-size must be >= 1.")
  }
  if (!is.finite(opt$ambiguous_eps_ratio) || opt$ambiguous_eps_ratio < 1) {
    stop("--ambiguous-eps-ratio must be >= 1.")
  }
  if (is.na(opt$plot_label_max_cluster_size) || opt$plot_label_max_cluster_size < 0) {
    stop("--plot-label-max-cluster-size must be >= 0.")
  }
  opt
}

output_path <- function(output_file, suffix) {
  sub("\\.json$", suffix, output_file)
}

parse_cli <- function(args) {
  if (length(args) == 0 || args[[1]] %in% c("--help", "-h")) {
    usage()
    quit(status = ifelse(length(args) >= 1 && args[[1]] %in% c("--help", "-h"), 0, 1))
  }

  first_option <- which(grepl("^--", args))[1]
  if (is.na(first_option)) {
    positional <- args
    extra_args <- character()
  } else {
    positional <- args[seq_len(first_option - 1)]
    extra_args <- args[first_option:length(args)]
  }

  if (!length(positional) %in% c(2, 3, 5)) {
    usage()
    stop(
      "Expected either input/output, input/output/similarity_threshold, ",
      "or the legacy five positional arguments."
    )
  }

  opt <- parse_options(extra_args)
  opt$input_file <- positional[[1]]
  opt$output_file <- positional[[2]]

  if (length(positional) >= 3) {
    opt$similarity_threshold <- positional[[3]]
  }
  if (length(positional) == 5) {
    opt$levels <- as.integer(positional[[5]])
  }
  if (is.na(opt$levels) || opt$levels < 1) {
    stop("levels must be a positive integer.")
  }

  opt
}

safe_mean <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  mean(values)
}

safe_min <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  min(values)
}

safe_max <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  max(values)
}

entropy_from_counts <- function(counts) {
  total <- sum(counts)
  if (total <= 0) {
    return(NA_real_)
  }
  probs <- counts[counts > 0] / total
  -sum(probs * log(probs))
}

read_distance_matrices <- function(input_file) {
  df <- fread(input_file, header = TRUE)
  regular_matrix <- acast(df, group.a ~ group.b, value.var = "estimated.difference.rate")
  regular_matrix <- as.matrix(regular_matrix)
  storage.mode(regular_matrix) <- "numeric"

  if (nrow(regular_matrix) != ncol(regular_matrix)) {
    stop("Input distance table did not produce a square group.a/group.b matrix.")
  }
  if (!all(rownames(regular_matrix) == colnames(regular_matrix))) {
    common <- intersect(rownames(regular_matrix), colnames(regular_matrix))
    regular_matrix <- regular_matrix[common, common, drop = FALSE]
  }
  diag(regular_matrix) <- 0

  finite_values <- regular_matrix[is.finite(regular_matrix)]
  finite_values <- finite_values[finite_values >= 0]
  if (length(finite_values) == 0) {
    stop("No finite distances found in input file.")
  }
  max_d <- max(finite_values)
  if (!is.finite(max_d) || max_d == 0) {
    max_d <- 1
  }

  norm_matrix <- regular_matrix / max_d
  diag(norm_matrix) <- 0
  if (any(!is.finite(norm_matrix))) {
    stop("Normalized distance matrix contains non-finite values.")
  }

  list(
    regular = regular_matrix,
    normalized = norm_matrix,
    max_distance = max_d,
    mpd_norm = safe_mean(norm_matrix[upper.tri(norm_matrix)])
  )
}

run_dbscan_eps <- function(distance_matrix, eps) {
  clustering <- dbscan(distance_matrix, eps = eps, minPts = 1)$cluster
  names(clustering) <- labels(distance_matrix)
  clustering
}

cluster_sizes <- function(clustering) {
  sort(table(as.character(clustering)), decreasing = TRUE)
}

pairwise_values_for_clusters <- function(matrix, clustering) {
  names_vec <- names(clustering)
  groups <- split(names_vec, as.character(clustering))
  group_names <- names(groups)
  within_values <- numeric()
  between_values <- numeric()

  for (group in group_names) {
    members <- groups[[group]]
    if (length(members) > 1) {
      sub <- matrix[members, members, drop = FALSE]
      within_values <- c(within_values, sub[upper.tri(sub)])
    }
  }

  if (length(group_names) > 1) {
    for (i in seq_len(length(group_names) - 1)) {
      for (j in (i + 1):length(group_names)) {
        members_i <- groups[[group_names[[i]]]]
        members_j <- groups[[group_names[[j]]]]
        between_values <- c(
          between_values,
          as.numeric(matrix[members_i, members_j, drop = FALSE])
        )
      }
    }
  }

  list(within = within_values, between = between_values)
}

cluster_nearest_distances <- function(matrix, clustering) {
  names_vec <- names(clustering)
  groups <- split(names_vec, as.character(clustering))
  group_names <- names(groups)

  rows <- lapply(group_names, function(group) {
    members <- groups[[group]]
    other_groups <- setdiff(group_names, group)

    nearest_cluster <- NA_character_
    nearest_distance_norm <- NA_real_
    if (length(other_groups) > 0) {
      distances <- vapply(other_groups, function(other) {
        other_members <- groups[[other]]
        mean(matrix[members, other_members, drop = FALSE])
      }, numeric(1))
      nearest_cluster <- names(distances)[which.min(distances)]
      nearest_distance_norm <- min(distances)
    }

    data.frame(
      cluster = group,
      cluster_size = length(members),
      nearest_cluster = nearest_cluster,
      nearest_distance_norm = nearest_distance_norm
    )
  })

  rbindlist(rows, fill = TRUE)
}

order_clustering_for_distance <- function(clustering, distance_matrix) {
  distance_labels <- labels(distance_matrix)
  if (is.null(distance_labels)) {
    return(clustering)
  }
  if (is.null(names(clustering))) {
    if (length(clustering) != length(distance_labels)) {
      stop("Unnamed clustering length does not match distance matrix size.")
    }
    names(clustering) <- distance_labels
    return(clustering)
  }
  missing <- setdiff(distance_labels, names(clustering))
  if (length(missing) > 0) {
    stop(
      "Clustering is missing label(s) present in the distance matrix: ",
      paste(missing, collapse = ", ")
    )
  }
  clustering[distance_labels]
}

compute_dbcv_score <- function(distance_matrix, clustering, dbcv_dim = 2) {
  n_clusters <- length(unique(clustering))
  if (length(clustering) < 3 || n_clusters < 2 || n_clusters >= length(clustering)) {
    return(NA_real_)
  }

  result <- tryCatch(
    suppressWarnings(
      dbscan::dbcv(distance_matrix, as.integer(factor(clustering)), d = dbcv_dim)
    ),
    error = function(e) NULL
  )
  if (is.null(result)) {
    return(NA_real_)
  }
  if (is.list(result) && "score" %in% names(result)) {
    score <- as.numeric(result$score)[[1]]
    return(ifelse(is.finite(score), score, NA_real_))
  }
  tryCatch(
    {
      score <- as.numeric(result)[[1]]
      ifelse(is.finite(score), score, NA_real_)
    },
    error = function(e) NA_real_
  )
}

compute_cluster_metrics <- function(clustering,
                                    norm_matrix,
                                    distance_matrix = as.dist(norm_matrix),
                                    eps = NA_real_,
                                    mpd_norm = NA_real_,
                                    compute_dbcv = FALSE,
                                    dbcv_dim = 2,
                                    small_cluster_size = 2,
                                    ambiguous_eps_ratio = 1.25) {
  clustering <- order_clustering_for_distance(clustering, distance_matrix)
  clustering <- setNames(as.character(clustering), names(clustering))
  n_haplotypes <- length(clustering)
  sizes <- as.integer(cluster_sizes(clustering))
  n_clusters <- length(sizes)
  singleton_count <- sum(sizes == 1)
  small_cluster_count <- sum(sizes <= small_cluster_size)
  largest_size <- if (length(sizes)) max(sizes) else NA_integer_
  pairwise <- pairwise_values_for_clusters(norm_matrix, clustering)
  nearest <- cluster_nearest_distances(norm_matrix, clustering)
  nearest$nearest_eps_ratio <- if (is.finite(eps) && eps > 0) {
    nearest$nearest_distance_norm / eps
  } else {
    NA_real_
  }
  small_nearest <- nearest[nearest$cluster_size <= small_cluster_size, ]
  singleton_nearest <- nearest[nearest$cluster_size == 1, ]
  ambiguous_small_cluster_count <- if (is.finite(eps) && eps > 0) {
    sum(
      small_nearest$nearest_eps_ratio <= ambiguous_eps_ratio,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }

  max_within <- safe_max(pairwise$within)
  min_between <- safe_min(pairwise$between)
  separation_ratio <- if (is.finite(max_within) && max_within > 0) {
    min_between / max_within
  } else {
    NA_real_
  }

  silhouette_mean <- NA_real_
  silhouette_min <- NA_real_
  if (n_haplotypes > 2 && n_clusters > 1 && n_clusters < n_haplotypes) {
    sil <- tryCatch(
      cluster::silhouette(as.integer(factor(clustering)), distance_matrix),
      error = function(e) NULL
    )
    if (!is.null(sil)) {
      silhouette_mean <- mean(sil[, "sil_width"])
      silhouette_min <- min(sil[, "sil_width"])
    }
  }

  dbcv_score <- if (compute_dbcv) {
    compute_dbcv_score(distance_matrix, clustering, dbcv_dim = dbcv_dim)
  } else {
    NA_real_
  }

  data.frame(
    eps = eps,
    num_haplotypes = n_haplotypes,
    num_clusters = n_clusters,
    singleton_cluster_count = singleton_count,
    singleton_cluster_fraction = singleton_count / n_clusters,
    small_cluster_count = small_cluster_count,
    small_cluster_fraction = small_cluster_count / n_clusters,
    ambiguous_small_cluster_count = ambiguous_small_cluster_count,
    ambiguous_small_cluster_fraction = ambiguous_small_cluster_count / n_clusters,
    largest_cluster_size = largest_size,
    largest_cluster_fraction = largest_size / n_haplotypes,
    cluster_size_entropy = entropy_from_counts(sizes),
    singleton_nearest_distance_norm_mean = safe_mean(singleton_nearest$nearest_distance_norm),
    singleton_nearest_distance_norm_min = safe_min(singleton_nearest$nearest_distance_norm),
    singleton_nearest_eps_ratio_mean = safe_mean(singleton_nearest$nearest_eps_ratio),
    singleton_nearest_eps_ratio_min = safe_min(singleton_nearest$nearest_eps_ratio),
    small_cluster_nearest_distance_norm_mean = safe_mean(small_nearest$nearest_distance_norm),
    small_cluster_nearest_distance_norm_min = safe_min(small_nearest$nearest_distance_norm),
    small_cluster_nearest_eps_ratio_mean = safe_mean(small_nearest$nearest_eps_ratio),
    small_cluster_nearest_eps_ratio_min = safe_min(small_nearest$nearest_eps_ratio),
    mean_within_distance_norm = safe_mean(pairwise$within),
    max_within_distance_norm = max_within,
    mean_between_distance_norm = safe_mean(pairwise$between),
    min_between_distance_norm = min_between,
    separation_ratio = separation_ratio,
    silhouette_mean = silhouette_mean,
    silhouette_min = silhouette_min,
    dbcv_score = dbcv_score,
    mpd_norm = mpd_norm,
    region_similarity = 1 - mpd_norm
  )
}

bounded_component <- function(value, lower = -1, upper = 1) {
  if (!is.finite(value)) {
    return(0)
  }
  min(upper, max(lower, value))
}

quality_score <- function(metrics_row, opt) {
  n_clusters <- metrics_row$num_clusters
  n_haplotypes <- metrics_row$num_haplotypes
  mpd <- metrics_row$mpd_norm
  low_diversity <- is.finite(mpd) && mpd <= opt$low_diversity_mpd_norm

  if (n_clusters == 1) {
    if (low_diversity) {
      if (opt$low_diversity_mpd_norm == 0) {
        return(0.75)
      }
      return(0.75 + 0.25 * (opt$low_diversity_mpd_norm - mpd) /
        opt$low_diversity_mpd_norm)
    }
    return(-0.75 - bounded_component(mpd, lower = 0, upper = 1))
  }

  silhouette_component <- bounded_component(metrics_row$silhouette_mean)
  dbcv_component <- bounded_component(metrics_row$dbcv_score)
  separation_component <- if (is.finite(metrics_row$separation_ratio)) {
    min(log1p(metrics_row$separation_ratio) / log1p(10), 1)
  } else {
    0
  }

  score <- silhouette_component + 0.50 * dbcv_component + 0.15 * separation_component

  singleton_fraction <- metrics_row$singleton_cluster_fraction
  ambiguous_small_fraction <- metrics_row$ambiguous_small_cluster_fraction
  if (is.finite(singleton_fraction)) {
    score <- score - 0.08 * singleton_fraction
  }
  if (is.finite(ambiguous_small_fraction)) {
    score <- score - 0.20 * ambiguous_small_fraction
  }

  # Singleton/small clusters are acceptable if they are clearly isolated.
  # Reward small clusters whose nearest other cluster is well beyond eps.
  isolation_ratio <- metrics_row$small_cluster_nearest_eps_ratio_mean
  if (is.finite(isolation_ratio) && is.finite(singleton_fraction)) {
    isolation_excess <- max(0, isolation_ratio - opt$ambiguous_eps_ratio)
    isolation_credit <- min(isolation_excess / opt$ambiguous_eps_ratio, 1)
    score <- score + 0.06 * isolation_credit * singleton_fraction
  }

  largest_fraction <- metrics_row$largest_cluster_fraction
  if (!low_diversity && is.finite(largest_fraction) &&
      largest_fraction > opt$giant_cluster_fraction) {
    denominator <- 1 - opt$giant_cluster_fraction
    giant_excess <- (largest_fraction - opt$giant_cluster_fraction) / denominator
    score <- score - 0.35 * giant_excess
  }

  if (n_clusters == n_haplotypes) {
    score <- score - 0.25
  }

  score
}

candidate_eps_values <- function(opt) {
  values <- seq(opt$eps_min, opt$eps_max, by = opt$eps_step)
  values <- unique(c(values, opt$eps_max))
  sort(values)
}

find_legacy_eps <- function(distance_matrix, region_similarity, opt) {
  optimal_eps <- opt$eps_max
  eps_values <- candidate_eps_values(opt)
  previous_clusters <- length(table(run_dbscan_eps(distance_matrix, eps_values[[1]])))

  for (eps in eps_values[-1]) {
    current_clusters <- length(table(run_dbscan_eps(distance_matrix, eps)))
    if (abs(previous_clusters - current_clusters) <= 1) {
      if ((current_clusters <= round(attr(distance_matrix, "Size") / 10)) ||
          region_similarity < 0.9) {
        optimal_eps <- eps
        break
      }
    }
    previous_clusters <- current_clusters
  }
  optimal_eps
}

choose_eps <- function(distance_matrix,
                       similarity_threshold,
                       norm_matrix,
                       mpd_norm,
                       opt) {
  if (similarity_threshold != "automatic") {
    eps <- 1 - as.numeric(similarity_threshold)
    if (!is.finite(eps)) {
      stop("similarity_threshold must be 'automatic' or numeric.")
    }
    return(list(
      eps = max(0, eps),
      method = "fixed_similarity_threshold",
      scan = data.frame()
    ))
  }

  if (opt$eps_selection == "legacy") {
    eps <- find_legacy_eps(distance_matrix, 1 - mpd_norm, opt)
    return(list(
      eps = eps,
      method = "legacy_stability_scan",
      scan = data.frame()
    ))
  }

  scan_rows <- lapply(candidate_eps_values(opt), function(eps) {
    clustering <- run_dbscan_eps(distance_matrix, eps)
    metrics <- compute_cluster_metrics(
      clustering,
      norm_matrix,
      distance_matrix = distance_matrix,
      eps = eps,
      mpd_norm = mpd_norm,
      compute_dbcv = opt$score_use_dbcv,
      dbcv_dim = opt$dbcv_dim,
      small_cluster_size = opt$small_cluster_size,
      ambiguous_eps_ratio = opt$ambiguous_eps_ratio
    )
    metrics$quality_score <- quality_score(metrics, opt)
    metrics
  })
  scan <- rbindlist(scan_rows, fill = TRUE)
  scan <- scan[order(
    -quality_score,
    -silhouette_mean,
    -dbcv_score,
    -separation_ratio,
    ambiguous_small_cluster_fraction,
    singleton_cluster_fraction,
    largest_cluster_fraction,
    num_clusters
  )]
  eps <- scan$eps[[1]]
  scan <- scan[order(eps)]
  scan$selected <- scan$eps == eps

  list(
    eps = eps,
    method = "quality_scan",
    scan = as.data.frame(scan)
  )
}

cluster_recursive <- function(names_vec,
                              norm_dist_mat,
                              regular_mat,
                              level,
                              prefix,
                              max_level,
                              results,
                              parent_eps,
                              opt) {
  norm_mat <- as.matrix(norm_dist_mat)
  local_mpd <- safe_mean(norm_mat[upper.tri(norm_mat)])

  if (level == 1) {
    current_eps <- parent_eps
  } else {
    local_choice <- choose_eps(
      norm_dist_mat,
      "automatic",
      norm_mat,
      local_mpd,
      opt
    )
    current_eps <- local_choice$eps
    if (!is.finite(current_eps) || current_eps <= 0 || current_eps > opt$eps_max) {
      current_eps <- max(parent_eps * 0.75, 0.05)
    }
  }

  clustering <- run_dbscan_eps(norm_dist_mat, current_eps)
  names(clustering) <- names_vec

  for (cl in sort(unique(clustering))) {
    members <- names(clustering[clustering == cl])
    group_name <- ifelse(prefix == "", paste0("HaploGroup", cl), paste0(prefix, ".", cl))

    for (m in members) {
      results[[m]] <- group_name
    }

    if (level < max_level && length(members) > 2) {
      raw_submatrix <- regular_mat[members, members, drop = FALSE]
      local_max <- max(raw_submatrix[is.finite(raw_submatrix)])

      if (local_max == 0 || is.infinite(local_max)) {
        next
      }

      local_norm <- raw_submatrix / local_max
      diag(local_norm) <- 0
      sub_dist <- as.dist(local_norm)

      results <- cluster_recursive(
        members,
        sub_dist,
        raw_submatrix,
        level + 1,
        group_name,
        max_level,
        results,
        current_eps,
        opt
      )
    }
  }

  results
}

cluster_distance_matrices <- function(final_vec, norm_matrix, regular_matrix) {
  group_names <- unique(unname(final_vec))
  k <- length(group_names)

  cluster_dist_norm <- matrix(0, nrow = k, ncol = k)
  rownames(cluster_dist_norm) <- group_names
  colnames(cluster_dist_norm) <- group_names
  cluster_dist <- cluster_dist_norm

  if (k > 1) {
    for (i in seq_len(k - 1)) {
      for (j in (i + 1):k) {
        g1 <- group_names[[i]]
        g2 <- group_names[[j]]
        members_i <- names(final_vec[final_vec == g1])
        members_j <- names(final_vec[final_vec == g2])
        if (length(members_i) > 0 && length(members_j) > 0) {
          dnorm <- norm_matrix[members_i, members_j, drop = FALSE]
          dreg <- regular_matrix[members_i, members_j, drop = FALSE]
          cluster_dist_norm[i, j] <- cluster_dist_norm[j, i] <- mean(dnorm)
          cluster_dist[i, j] <- cluster_dist[j, i] <- mean(dreg)
        }
      }
    }
  }

  list(norm = cluster_dist_norm, raw = cluster_dist)
}

cluster_summary <- function(final_vec, norm_matrix, regular_matrix) {
  groups <- split(names(final_vec), unname(final_vec))
  group_names <- names(groups)

  rows <- lapply(group_names, function(group) {
    members <- groups[[group]]

    mean_within_norm <- NA_real_
    max_within_norm <- NA_real_
    mean_within_raw <- NA_real_
    max_within_raw <- NA_real_
    medoid <- members[[1]]

    if (length(members) > 1) {
      dnorm <- norm_matrix[members, members, drop = FALSE]
      draw <- regular_matrix[members, members, drop = FALSE]
      mean_norm_by_member <- rowMeans(dnorm)
      medoid <- names(mean_norm_by_member[which.min(mean_norm_by_member)])
      mean_within_norm <- safe_mean(dnorm[upper.tri(dnorm)])
      max_within_norm <- safe_max(dnorm[upper.tri(dnorm)])
      mean_within_raw <- safe_mean(draw[upper.tri(draw)])
      max_within_raw <- safe_max(draw[upper.tri(draw)])
    }

    other_groups <- setdiff(group_names, group)
    nearest_cluster <- NA_character_
    nearest_cluster_distance_norm <- NA_real_
    nearest_cluster_distance_raw <- NA_real_
    if (length(other_groups) > 0) {
      distances <- lapply(other_groups, function(other) {
        other_members <- groups[[other]]
        data.frame(
          nearest_cluster = other,
          nearest_cluster_distance_norm = mean(norm_matrix[members, other_members, drop = FALSE]),
          nearest_cluster_distance_raw = mean(regular_matrix[members, other_members, drop = FALSE])
        )
      })
      distances <- rbindlist(distances)
      distances <- distances[order(nearest_cluster_distance_norm)]
      nearest_cluster <- distances$nearest_cluster[[1]]
      nearest_cluster_distance_norm <- distances$nearest_cluster_distance_norm[[1]]
      nearest_cluster_distance_raw <- distances$nearest_cluster_distance_raw[[1]]
    }

    data.frame(
      haplotype.group = group,
      cluster_size = length(members),
      medoid = medoid,
      mean_within_distance_norm = mean_within_norm,
      max_within_distance_norm = max_within_norm,
      mean_within_distance = mean_within_raw,
      max_within_distance = max_within_raw,
      nearest_cluster = nearest_cluster,
      nearest_cluster_distance_norm = nearest_cluster_distance_norm,
      nearest_cluster_distance = nearest_cluster_distance_raw,
      haplotypes = paste(members, collapse = ",")
    )
  })

  rbindlist(rows, fill = TRUE)
}

haplotype_cluster_diagnostics <- function(final_vec, norm_matrix, regular_matrix, eps) {
  groups <- split(names(final_vec), unname(final_vec))
  group_names <- names(groups)

  rows <- lapply(names(final_vec), function(haplotype) {
    group <- unname(final_vec[[haplotype]])
    members <- groups[[group]]
    own_members <- setdiff(members, haplotype)

    own_mean_norm <- NA_real_
    own_min_norm <- NA_real_
    own_mean_raw <- NA_real_
    own_min_raw <- NA_real_
    if (length(own_members) > 0) {
      own_mean_norm <- mean(norm_matrix[haplotype, own_members])
      own_min_norm <- min(norm_matrix[haplotype, own_members])
      own_mean_raw <- mean(regular_matrix[haplotype, own_members])
      own_min_raw <- min(regular_matrix[haplotype, own_members])
    }

    other_groups <- setdiff(group_names, group)
    nearest_cluster <- NA_character_
    nearest_cluster_mean_distance_norm <- NA_real_
    nearest_cluster_min_distance_norm <- NA_real_
    nearest_cluster_mean_distance <- NA_real_
    nearest_cluster_min_distance <- NA_real_
    nearest_haplotype <- NA_character_

    if (length(other_groups) > 0) {
      distances <- lapply(other_groups, function(other) {
        other_members <- groups[[other]]
        norm_values <- norm_matrix[haplotype, other_members]
        raw_values <- regular_matrix[haplotype, other_members]
        data.frame(
          nearest_cluster = other,
          mean_distance_norm = mean(norm_values),
          min_distance_norm = min(norm_values),
          mean_distance = mean(raw_values),
          min_distance = min(raw_values),
          nearest_haplotype = other_members[[which.min(norm_values)]]
        )
      })
      distances <- rbindlist(distances)
      distances <- distances[order(mean_distance_norm, min_distance_norm)]
      nearest_cluster <- distances$nearest_cluster[[1]]
      nearest_cluster_mean_distance_norm <- distances$mean_distance_norm[[1]]
      nearest_cluster_min_distance_norm <- distances$min_distance_norm[[1]]
      nearest_cluster_mean_distance <- distances$mean_distance[[1]]
      nearest_cluster_min_distance <- distances$min_distance[[1]]
      nearest_haplotype <- distances$nearest_haplotype[[1]]
    }

    nearest_eps_ratio <- if (is.finite(eps) && eps > 0) {
      nearest_cluster_mean_distance_norm / eps
    } else {
      NA_real_
    }

    data.frame(
      haplotype.name = haplotype,
      haplotype.group = group,
      cluster_size = length(members),
      mean_distance_to_assigned_cluster_norm = own_mean_norm,
      min_distance_to_assigned_cluster_norm = own_min_norm,
      mean_distance_to_assigned_cluster = own_mean_raw,
      min_distance_to_assigned_cluster = own_min_raw,
      nearest_other_cluster = nearest_cluster,
      nearest_other_cluster_mean_distance_norm = nearest_cluster_mean_distance_norm,
      nearest_other_cluster_min_distance_norm = nearest_cluster_min_distance_norm,
      nearest_other_cluster_mean_distance = nearest_cluster_mean_distance,
      nearest_other_cluster_min_distance = nearest_cluster_min_distance,
      nearest_other_haplotype = nearest_haplotype,
      nearest_other_cluster_eps_ratio = nearest_eps_ratio
    )
  })

  rbindlist(rows, fill = TRUE)
}

write_cluster_plot <- function(output_file,
                               final_vec,
                               norm_matrix,
                               opt,
                               eps,
                               mpd_norm) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not available; skipping cluster plot.")
    return(invisible(FALSE))
  }
  if (length(final_vec) < 2) {
    warning("Fewer than two haplotypes; skipping cluster plot.")
    return(invisible(FALSE))
  }

  distance_matrix <- as.dist(norm_matrix)
  pcoa <- tryCatch(
    cmdscale(distance_matrix, k = 2, eig = TRUE),
    error = function(e) NULL
  )
  if (is.null(pcoa)) {
    warning("Could not compute PCoA coordinates; skipping cluster plot.")
    return(invisible(FALSE))
  }

  points <- as.data.frame(pcoa$points)
  if (ncol(points) == 1) {
    points$V2 <- 0
  }
  colnames(points)[1:2] <- c("pc1", "pc2")
  points$haplotype <- rownames(points)
  points$haplotype.group <- unname(final_vec[points$haplotype])

  sizes <- as.data.frame(table(unname(final_vec)))
  colnames(sizes) <- c("haplotype.group", "cluster_size")
  points <- merge(points, sizes, by = "haplotype.group", all.x = TRUE)
  points$label <- ifelse(
    points$cluster_size <= opt$plot_label_max_cluster_size,
    points$haplotype,
    NA_character_
  )

  title <- paste0(
    "COSIGT haplogroup clusters (eps=", signif(eps, 4),
    ", mpd_norm=", signif(mpd_norm, 4), ")"
  )
  plot <- ggplot2::ggplot(
    points,
    ggplot2::aes(x = pc1, y = pc2, color = haplotype.group)
  ) +
    ggplot2::geom_point(ggplot2::aes(size = cluster_size), alpha = 0.85) +
    ggplot2::scale_size_continuous(range = c(1.8, 5)) +
    ggplot2::labs(
      x = "PCoA1",
      y = "PCoA2",
      color = "Cluster",
      size = "Cluster size",
      title = title
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  label_data <- points[!is.na(points$label), ]
  if (nrow(label_data) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = label),
        size = 2.5,
        max.overlaps = Inf,
        show.legend = FALSE
      )
    } else {
      plot <- plot + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = label),
        size = 2.5,
        vjust = -0.6,
        show.legend = FALSE
      )
    }
  }

  png_path <- output_path(output_file, ".cluster_pcoa.png")
  pdf_path <- output_path(output_file, ".cluster_pcoa.pdf")
  ggplot2::ggsave(
    png_path,
    plot,
    width = opt$plot_width,
    height = opt$plot_height,
    dpi = 300,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = opt$plot_width,
    height = opt$plot_height,
    limitsize = FALSE
  )
  invisible(TRUE)
}

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
input_file <- opt$input_file
output_file <- opt$output_file
similarity_threshold <- opt$similarity_threshold
levels <- opt$levels

matrices <- read_distance_matrices(input_file)
regular_matrix <- matrices$regular
norm_matrix <- matrices$normalized
mpd_norm <- matrices$mpd_norm
distance_matrix <- as.dist(norm_matrix)

eps_choice <- choose_eps(
  distance_matrix,
  similarity_threshold,
  norm_matrix,
  mpd_norm,
  opt
)
eps <- eps_choice$eps

final_res <- cluster_recursive(
  names_vec = labels(distance_matrix),
  norm_dist_mat = distance_matrix,
  regular_mat = regular_matrix,
  level = 1,
  prefix = "",
  max_level = levels,
  results = list(),
  parent_eps = eps,
  opt = opt
)
final_vec <- unlist(final_res)

write(toJSON(final_res), output_file)

reversed_data <- split(names(final_vec), unname(final_vec))
haplotable <- data.frame(
  haplotype.name = unlist(reversed_data),
  haplotype.group = rep(names(reversed_data), lengths(reversed_data))
)
rownames(haplotable) <- NULL
fwrite(
  haplotable,
  output_path(output_file, ".tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

cluster_distances <- cluster_distance_matrices(final_vec, norm_matrix, regular_matrix)
fwrite(
  data.frame(h.group = row.names(cluster_distances$norm), cluster_distances$norm),
  output_path(output_file, ".hapdist.norm.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)
fwrite(
  data.frame(h.group = row.names(cluster_distances$raw), cluster_distances$raw),
  output_path(output_file, ".hapdist.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

summary_table <- cluster_summary(final_vec, norm_matrix, regular_matrix)
fwrite(
  summary_table,
  output_path(output_file, ".cluster_summary.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

haplotype_diagnostics <- haplotype_cluster_diagnostics(
  final_vec,
  norm_matrix,
  regular_matrix,
  eps
)
fwrite(
  haplotype_diagnostics,
  output_path(output_file, ".haplotype_cluster_diagnostics.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

medoids <- summary_table[, .(haplotype.group, medoid)]
fwrite(
  medoids,
  output_path(output_file, ".medoids.tsv"),
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)

metrics <- compute_cluster_metrics(
  final_vec,
  norm_matrix,
  distance_matrix = distance_matrix,
  eps = eps,
  mpd_norm = mpd_norm,
  compute_dbcv = TRUE,
  dbcv_dim = opt$dbcv_dim,
  small_cluster_size = opt$small_cluster_size,
  ambiguous_eps_ratio = opt$ambiguous_eps_ratio
)
metrics$method <- "dbscan_minpts1_all_assigned"
metrics$eps_selection <- eps_choice$method
metrics$max_raw_distance <- matrices$max_distance
metrics$levels <- levels
metrics$score_use_dbcv <- opt$score_use_dbcv
metrics$dbcv_dim <- opt$dbcv_dim
metrics$low_diversity_mpd_norm <- opt$low_diversity_mpd_norm
metrics$giant_cluster_fraction_threshold <- opt$giant_cluster_fraction
metrics$small_cluster_size_threshold <- opt$small_cluster_size
metrics$ambiguous_eps_ratio_threshold <- opt$ambiguous_eps_ratio
metrics$quality_score <- quality_score(metrics, opt)
metrics$eps_minus_0.01_num_clusters <- length(table(run_dbscan_eps(
  distance_matrix,
  max(opt$eps_min, eps - 0.01)
)))
metrics$eps_plus_0.01_num_clusters <- length(table(run_dbscan_eps(
  distance_matrix,
  min(opt$eps_max, eps + 0.01)
)))
metrics$eps_minus_0.02_num_clusters <- length(table(run_dbscan_eps(
  distance_matrix,
  max(opt$eps_min, eps - 0.02)
)))
metrics$eps_plus_0.02_num_clusters <- length(table(run_dbscan_eps(
  distance_matrix,
  min(opt$eps_max, eps + 0.02)
)))

metric_order <- c(
  "method",
  "eps_selection",
  "eps",
  "levels",
  "num_haplotypes",
  "num_clusters",
  "mpd_norm",
  "region_similarity",
  "max_raw_distance",
  "score_use_dbcv",
  "dbcv_dim",
  "low_diversity_mpd_norm",
  "giant_cluster_fraction_threshold",
  "small_cluster_size_threshold",
  "ambiguous_eps_ratio_threshold",
  "quality_score",
  "singleton_cluster_count",
  "singleton_cluster_fraction",
  "small_cluster_count",
  "small_cluster_fraction",
  "ambiguous_small_cluster_count",
  "ambiguous_small_cluster_fraction",
  "largest_cluster_size",
  "largest_cluster_fraction",
  "cluster_size_entropy",
  "singleton_nearest_distance_norm_mean",
  "singleton_nearest_distance_norm_min",
  "singleton_nearest_eps_ratio_mean",
  "singleton_nearest_eps_ratio_min",
  "small_cluster_nearest_distance_norm_mean",
  "small_cluster_nearest_distance_norm_min",
  "small_cluster_nearest_eps_ratio_mean",
  "small_cluster_nearest_eps_ratio_min",
  "mean_within_distance_norm",
  "max_within_distance_norm",
  "mean_between_distance_norm",
  "min_between_distance_norm",
  "separation_ratio",
  "silhouette_mean",
  "silhouette_min",
  "dbcv_score",
  "eps_minus_0.01_num_clusters",
  "eps_plus_0.01_num_clusters",
  "eps_minus_0.02_num_clusters",
  "eps_plus_0.02_num_clusters"
)
metrics <- metrics[, metric_order]
fwrite(
  metrics,
  output_path(output_file, ".metrics.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

if (nrow(eps_choice$scan) > 0) {
  fwrite(
    eps_choice$scan,
    output_path(output_file, ".eps_scan.tsv"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
}

if (isTRUE(opt$plot_clusters)) {
  write_cluster_plot(output_file, final_vec, norm_matrix, opt, eps, mpd_norm)
}
