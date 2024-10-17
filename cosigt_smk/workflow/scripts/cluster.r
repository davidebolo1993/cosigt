#!/usr/bin/env Rscript 

# Load required libraries
library(data.table)
library(reshape2)
library(NbClust)
library(rjson)
library(dendextend)
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Set data.table threads
setDTthreads(1)

# Read and process input data
df <- fread(input_file)
df$jaccard.distance <- 1 - df$jaccard.similarity

# Create distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "jaccard.distance")
distanceMatrix <- as.dist(regularMatrix)

# Calculate silhouette score and best partition
max_clusters <- round(length(unique(df$group.a)) / 3)
res <- NbClust(diss = distanceMatrix, method = "average", index = "silhouette", 
               distance = NULL, max.nc = max_clusters)$Best.partition

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
distanceMatrix <- as.dist(regularMatrix)

# Perform hierarchical clustering
k <- as.numeric(tail(sort(res), 1))
hc <- hclust(distanceMatrix, method = "average")
dend <- as.dendrogram(hc) %>%
  set("branches_k_color", k = k) %>%
  set("labels_cex", 0.6)

# Create and save dendrogram plot
ggd1 <- as.ggdend(dend)
p <- ggplot(ggd1, horiz = TRUE) + theme_minimal()
ggsave(gsub(".json", ".pdf", output_file), plot = p, width = 20, height = 15)