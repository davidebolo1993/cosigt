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
region_similarity <- round(as.numeric(args[4]), 2) #ignored at the moment
levels <- as.integer(args[5])

df <- fread(input_file, header=TRUE)

# Distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "estimated.difference.rate")
maxD <- max(regularMatrix[is.finite(regularMatrix)])
normRegularMatrix <- regularMatrix / maxD#sort of global similarity
mpd_norm <- mean(normRegularMatrix[upper.tri(normRegularMatrix)])
region_similarity<-1-mpd_norm

# Function to determine optimal eps
find_optimal_eps <- function(distanceMatrix, similarity_threshold) {
  if (similarity_threshold != "automatic") {
    return(1 - as.numeric(similarity_threshold))
  }
  optimal_eps <- 0
  pclust <- length(table(dbscan(distanceMatrix, eps = 0, minPts = 1)$cluster))
  for (eps in seq(0.01, 0.30, 0.01)) {
    cclust <- length(table(dbscan(distanceMatrix, eps=eps, minPts=1)$cluster))
    if (abs(pclust - cclust) <= 1) {
      if ((cclust <= round(attr(distanceMatrix, "Size") / 10)) || region_similarity < 0.9) {
        optimal_eps <- eps
        break
      }
    }
    pclust <- cclust
  }
  return(ifelse(eps < 0.3, optimal_eps, 0.3))
}

cluster_recursive <- function(names_vec, norm_dist_mat, regular_mat, level, prefix, max_level, results, parent_eps) {
  if (level == 1) {
    current_eps <- parent_eps
  } else {
    # Subclustering: try to find optimal eps for this subgroup
    local_eps <- find_optimal_eps(norm_dist_mat, "automatic")
    if (local_eps == 0 || local_eps > 0.3) {
      current_eps <- parent_eps * 0.75
      current_eps <- max(current_eps, 0.05)
    } else {
      current_eps <- local_eps
    }
  }

  clustering <- dbscan(norm_dist_mat, eps=current_eps, minPts=1)$cluster
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

      # Recursively call with current_eps as parent_eps for next level
      results <- cluster_recursive(members, sub_dist, raw_submatrix, 
                                   level + 1, group_name, max_level, 
                                   results, current_eps)
    }
  }

  return(results)
}

# Top-level clustering
distanceMatrix <- as.dist(normRegularMatrix)
eps <- find_optimal_eps(distanceMatrix, similarity_threshold)

final_res <- cluster_recursive(names_vec=labels(distanceMatrix), 
                               norm_dist_mat=distanceMatrix, 
                               regular_mat=regularMatrix, 
                               level=1, 
                               prefix="", 
                               max_level=levels, 
                               results=list(), 
                               parent_eps=eps)

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
fwrite(haplotable, tsv_output, row.names = FALSE, col.names = TRUE, sep = "	")

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
fwrite(metrics, metrics_output, row.names=FALSE, col.names=TRUE, sep = "	")

distance_output <- gsub(".json", ".hapdist.norm.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist_norm), cluster_dist_norm), distance_output, row.names = FALSE, col.names = TRUE, sep = "	")

distance_output <- gsub(".json", ".hapdist.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist), cluster_dist), distance_output, row.names = FALSE, col.names = TRUE, sep = "	")

# Identify representative haplotype (medoid) for each cluster
center_hap <- matrix(0, nrow=k, ncol=2)
for (i in 1:k) {
  g1 <- group_names[i]
  members_i <- names(final_res[final_res == g1])
  if (length(members_i) > 0) {
    dnorm <- normRegularMatrix[members_i, members_i, drop = FALSE]
    mean_norm <- rowMeans(dnorm)
    center_hap[i,] <- c(g1, names(mean_norm[which.min(mean_norm)]))
  }
}

medoids <- gsub(".json", ".medoids.tsv", output_file)
fwrite(center_hap, medoids, row.names = FALSE, col.names = FALSE, sep = "	")
