#!/usr/bin/env Rscript 
library(data.table)
library(reshape2)
library(rjson)
library(dbscan)
library(cluster)

setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
similarity_threshold <- args[3]
region_similarity <- round(as.numeric(args[4]), 2)
levels <- ifelse(length(args) >= 5, as.integer(args[5]), 1)

df <- fread(input_file, header=TRUE)

# Distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "estimated.difference.rate")
maxD <- max(regularMatrix[is.finite(regularMatrix)])
normRegularMatrix <- regularMatrix / maxD

# Function to determine optimal eps
find_optimal_eps <- function(distanceMatrix, region_similarity, similarity_threshold) {
  if (similarity_threshold != "automatic") {
    return(1 - as.numeric(similarity_threshold))
  }
  optimal_eps <- 0
  pclust <- length(table(dbscan(distanceMatrix, eps = 0, minPts = 1)$cluster))
  for (eps in seq(0.01, 0.30, 0.01)) {
    cclust <- length(table(dbscan(distanceMatrix, eps=eps, minPts=1)$cluster))
    if (abs(pclust - cclust) <= 1) {
      if ((region_similarity >= 0.9 && cclust <= round(attr(distanceMatrix, "Size") / 10)) || 
          region_similarity < 0.9) {
        optimal_eps <- eps
        break
      }
    }
    pclust <- cclust
  }
  return(optimal_eps)
}

# Recursive clustering - the idea here is to have some
# sort of more accurate clustering 
cluster_recursive <- function(names_vec, norm_dist_mat, regular_mat, level, prefix, max_level, results, eps) {
  clustering <- dbscan(norm_dist_mat, eps=eps, minPts=1)$cluster
  names(clustering) <- names_vec
  for (cl in sort(unique(clustering))) {
    members <- names(clustering[clustering == cl])
    group_name <- ifelse(prefix == "", paste0("HaploGroup", cl), paste0(prefix, ".", cl))
    for (m in members) {
      results[[m]] <- group_name
    }
    if (level < max_level && length(members) > 2) {
      raw_submatrix <- regular_mat[members, members]
      local_max <- max(raw_submatrix[is.finite(raw_submatrix)])
      if (local_max == 0 || is.infinite(local_max)) {
        next
      }
      local_norm <- raw_submatrix / local_max
      sub_dist <- as.dist(local_norm)      
      sub_eps <- 0.2
      results <- cluster_recursive(members, sub_dist, raw_submatrix, level + 1, group_name, max_level, results, sub_eps)
    }
  }
  return(results)
}

# Top-level clustering
distanceMatrix <- as.dist(normRegularMatrix)
eps <- find_optimal_eps(distanceMatrix, region_similarity, similarity_threshold)
final_res <- cluster_recursive(labels(distanceMatrix), distanceMatrix, regularMatrix, level=1, prefix="", max_level=levels, results=list(), eps=eps)

# Write JSON
jout <- toJSON(final_res)
write(jout, output_file)

# Invert mapping
reversed_data <- list()
for (key in names(final_res)) {
  value <- final_res[[key]]
  reversed_data[[value]] <- c(reversed_data[[value]], key)
}

# Create haplotype table
haplotable <- data.frame(
  haplotype.name = unlist(reversed_data),
  haplotype.group = rep(names(reversed_data), lengths(reversed_data))
)
rownames(haplotable) <- NULL
tsv_output <- gsub(".json", ".tsv", output_file)
fwrite(haplotable, tsv_output, row.names = FALSE, col.names = TRUE, sep = "\t")

# Compute distances
all_groups <- unique(unlist(final_res))
group_names <- unique(all_groups)
k <- length(group_names)

cluster_dist_norm <- matrix(0, nrow = k, ncol = k)
rownames(cluster_dist_norm) <- group_names
colnames(cluster_dist_norm) <- group_names
cluster_dist <- cluster_dist_norm

if (k > 1) {
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      g1 <- group_names[i]
      g2 <- group_names[j]
      members_i <- names(final_res[final_res == g1])
      members_j <- names(final_res[final_res == g2])
      if (length(members_i) > 0 && length(members_j) > 0) {
        dnorm <- normRegularMatrix[members_i, members_j, drop = FALSE]
        dreg <- regularMatrix[members_i, members_j, drop = FALSE]
        mean_norm <- mean(dnorm)
        mean_reg <- mean(dreg)
        cluster_dist_norm[i, j] <- cluster_dist_norm[j, i] <- mean_norm
        cluster_dist[i, j] <- cluster_dist[j, i] <- mean_reg
      }
    }
  }
}

# Write metrics and distances
metrics_output <- gsub(".json", ".metrics.tsv", output_file)
metrics <- data.frame(
  eps = eps,
  num_clusters = k
)
fwrite(metrics, metrics_output, row.names=FALSE, col.names=TRUE, sep = "\t")

distance_output <- gsub(".json", ".hapdist.norm.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist_norm), cluster_dist_norm), distance_output, row.names = FALSE, col.names = TRUE, sep = "\t")

distance_output <- gsub(".json", ".hapdist.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist), cluster_dist), distance_output, row.names = FALSE, col.names = TRUE, sep = "\t")
