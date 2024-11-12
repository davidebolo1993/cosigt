#!/usr/bin/env Rscript 

# Load required libraries
library(data.table)
library(reshape2)
library(NbClust)
library(rjson)
library(dendextend)
library(randomcoloR)

# Set data.table threads
setDTthreads(1)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read and process input data
df <- fread(input_file)
df$jaccard.distance <- 1 - df$jaccard.similarity

#find outliers
iqr<-IQR(df$group.a.length)
q1<-as.numeric(quantile(df$group.a.length,probs=c(0.25)))
q3<-as.numeric(quantile(df$group.a.length,probs=c(0.75)))
lb<-q1-iqr*3
ub<-q3+iqr*3
#keep only those in range
newdf<-df[group.a.length >= lb & group.a.length <= ub & group.b.length >= lb & group.b.length <= ub]
#track excluded
excluded<-unique(df$group.a[df$group.a%in%newdf$group.a == FALSE])
df<-newdf

# Create distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "jaccard.distance")
#max distance if NA
regularMatrix[is.na(regularMatrix)]<-1
distanceMatrix <- as.dist(regularMatrix)

# Calculate silhouette score and best partition
max_cluster <- round(length(unique(df$group.a)) / 3) ##control
res <- NbClust(diss = distanceMatrix, method = "average", index = "silhouette", 
               distance = NULL, max.nc = max_cluster)$Best.partition

# Format results
res.list <- lapply(split(res, names(res)), unname)
named_res <- lapply(res.list, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
for (e in excluded) {named_res[[e]] <- "*"}
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
palette <- distinctColorPalette(k)

dend <- as.dendrogram(hc) %>%
  set("branches_k_color", value=palette, k = k) %>%
  set("labels_col", value=palette,k = k) %>%
  set("labels_cex", 0.6)

#map leaf colors to haplogroups
leaf_colors <- get_leaves_branches_col(dend)
leaf_names <- get_leaves_attr(dend, "label")
cluster_map <- data.frame(
  leaf = as.character(leaf_names),
  cluster = clusters[leaf_names],
  color = leaf_colors
)
cluster_colors <- sapply(unique(clusters), function(cluster) {
  leaf_index <- which(clusters == cluster)[1]
  leaf_colors[leaf_index]
})
cluster_info <- data.frame(
  cluster = paste0("HaploGroup", unique(clusters)),
  color = cluster_colors
)

unique_clusters <- cluster_map[!duplicated(cluster_map$cluster), ]
#actual plot
pdf(gsub(".json", ".pdf", output_file), width = 20, height = 15)
dend %>% plot(horiz = TRUE)
legend("topleft",
       legend = paste0("HaploGroup", unique_clusters$cluster),
       fill = unique_clusters$color,
       cex = 0.6)
dev.off()