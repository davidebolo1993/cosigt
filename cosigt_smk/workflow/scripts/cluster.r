#!/usr/bin/env Rscript 

# Load required libraries
library(data.table)
library(reshape2)
library(NbClust)
library(rjson)
library(dendextend)
library(ggplot2)
library(ggdendro)

# Set data.table threads
setDTthreads(1)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read and process input data
df <- fread(input_file)
df$jaccard.distance <- 1 - df$jaccard.similarity

# Create distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "jaccard.distance")
#max distance if NA
regularMatrix[is.na(regularMatrix)]<-1
distanceMatrix <- as.dist(regularMatrix)


# Calculate silhouette score and best partition
max_cluster <- round(length(unique(df$group.a)) / 3)
res <- NbClust(diss = distanceMatrix, method = "average", index = "silhouette", 
               distance = NULL, max.nc = max_cluster)$Best.partition

# Format results
res.list <- lapply(split(res, names(res)), unname)
named_res <- lapply(res.list, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
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

# Simplify matrix names for hclust
simplify_names <- function(x) unlist(strsplit(x, ":"))[1]
colnames(regularMatrix) <- sapply(colnames(regularMatrix), simplify_names)
rownames(regularMatrix) <- sapply(rownames(regularMatrix), simplify_names)

# Perform hierarchical clustering
k <- as.numeric(tail(sort(res), 1))
hc <- hclust(as.dist(regularMatrix), method = "average")

#calculate average distances between clusters
clusters <- cutree(hc, k = k)

cluster_dist <- matrix(0, nrow = k, ncol = k)
rownames(cluster_dist) <- paste0("HaploGroup", 1:k)
colnames(cluster_dist) <- paste0("HaploGroup", 1:k)
# Calculate mean distances between clusters
for(i in 1:(k-1)) {
  for(j in (i+1):k) {
    # Get indices for each cluster
    cluster_i <- which(clusters == i)
    cluster_j <- which(clusters == j)  
    # Calculate mean distance between clusters
    distances <- regularMatrix[cluster_i, cluster_j, drop = FALSE]
    mean_dist <- mean(distances)
      
    # Store distances symmetrically
    cluster_dist[i,j] <- mean_dist
    cluster_dist[j,i] <- mean_dist
  }
}
distance_output <- gsub(".json", ".hapdist.tsv", output_file)
fwrite(data.frame(h.group=row.names(cluster_dist),cluster_dist), distance_output, row.names = FALSE, col.names = TRUE, sep = "\t")

#plot
dend <- as.dendrogram(hc) %>%
  set("branches_k_color", k = k) %>%
  set("labels_cex", 0.6)

# Convert dendrogram to ggdendro format for ggplot
ggd1 <- as.ggdend(dend)

# Extract leaf labels and create a mapping of each label to its HaploGroup
leaf_labels <- labels(dend)
cluster_mapping <- data.frame(label = leaf_labels, HaploGroup = paste0("HaploGroup", clusters[leaf_labels]))

# Plot with branches colored by HaploGroup and labels as haplotype names
p <- ggplot(ggd1$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = factor(cluster_mapping$HaploGroup[match(xend, 1:length(leaf_labels))])),
               size = 0.8) +
  geom_text(data = cluster_mapping, aes(x = match(label, leaf_labels), y = -0.05, label = label), 
            angle = 0, vjust = 0.5, hjust = 0.2, color = "black", size=3) +
  scale_color_manual(name = "HaploGroup", values = rainbow(k)) +
  labs(color = "HaploGroup") +
  theme(legend.position = "right")+
  theme_minimal() +
  guides(color = guide_legend(na.translate = FALSE))  +
  coord_flip()

# Save the plot with branch colors and haplotype labels
ggsave(gsub(".json", ".pdf", output_file), plot = p, width = 20, height = 15)
