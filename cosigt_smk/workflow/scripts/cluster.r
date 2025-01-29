#!/usr/bin/env Rscript 

# Load required libraries
library(data.table)
library(reshape2)
library(rjson)
library(dbscan)

# Set data.table threads
setDTthreads(1)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
similarity_threshold<-as.numeric(args[3])

# Read and process input data
df <- fread(input_file)

# Create distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "estimated.difference.rate")
maxD<-max(regularMatrix[!is.na(regularMatrix)])
regularMatrix[is.na(regularMatrix)]<-Inf
normRegularMatrix<-regularMatrix/maxD
distanceMatrix <- as.dist(normRegularMatrix)

eps<-1-similarity_threshold
res <- dbscan(distanceMatrix, eps=eps, minPts = 1)$cluster
names(res)<-labels(distanceMatrix)

# Format results
res.list <- lapply(split(res, names(res)), unname)
named_res <- lapply(res, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
jout <- toJSON(named_res)

# Write JSON output
write(jout, output_file)

# Create reversed data
reversed_data <- list()
for (key in names(named_res)) {
    value <- named_res[[key]]
    if (!is.null(reversed_data[[value]])) {
        reversed_data[[value]] <- c(reversed_data[[value]], key)
    } else {
        reversed_data[[value]] <- key
    }
}

# Create haplotype table
haplotable <- data.frame(
    haplotype.name = unlist(reversed_data),
    haplotype.group = rep(names(reversed_data), lengths(reversed_data))
)
rownames(haplotable) <- NULL

# Write haplotype table
tsv_output <- gsub(".json", ".tsv", output_file)
fwrite(haplotable, tsv_output, row.names = FALSE, col.names = TRUE, sep = "\t")

k <- as.numeric(tail(sort(res), 1))
#write distances?
cluster_dist_norm <- matrix(0, nrow = k, ncol = k)
rownames(cluster_dist_norm) <- paste0("HaploGroup", 1:k)
colnames(cluster_dist_norm) <- paste0("HaploGroup", 1:k)
cluster_dist<-cluster_dist_norm
# Calculate mean distances between clusters
for(i in 1:(k-1)) {
  for(j in (i+1):k) {
    # Get indices for each cluster
    cluster_i <- which(res == i)
    cluster_j <- which(res == j)  
    # Calculate mean distance between clusters
    distances <- normRegularMatrix[cluster_i, cluster_j, drop = FALSE]
    mean_dist <- mean(distances)
    # Store distances symmetrically
    cluster_dist_norm[i,j] <- mean_dist
    cluster_dist_norm[j,i] <- mean_dist
    distances <- regularMatrix[cluster_i, cluster_j, drop = FALSE]
    mean_dist <- mean(distances)
    cluster_dist[i,j] <- mean_dist
    cluster_dist[j,i] <- mean_dist
  }
}

distance_output <- gsub(".json", ".hapdist.norm.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist_norm),cluster_dist_norm), distance_output, row.names = FALSE, col.names = TRUE, sep = "\t")
distance_output <- gsub(".json", ".hapdist.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist),cluster_dist), distance_output, row.names = FALSE, col.names = TRUE, sep = "\t")