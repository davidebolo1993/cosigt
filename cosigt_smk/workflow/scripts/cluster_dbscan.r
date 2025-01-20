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

#find outliers
z.scores.a <- (df$group.a.length - mean(df$group.a.length)) / sd(df$group.a.length)
z.scores.b <- (df$group.b.length - mean(df$group.b.length)) / sd(df$group.b.length)
outliers.a<-df[abs(z.scores.a) > 3]
outliers.b<-df[abs(z.scores.b) > 3]
outliers<-unique(c(which(abs(z.scores.a) > 3),which(abs(z.scores.b) > 3)))
if (length(outliers)==0) {
  excluded<-c()
} else {
  newdf<-df[-outliers]
  excluded<-unique(df$group.a[df$group.a%in%newdf$group.a == FALSE])
  df<-newdf
}

# Create distance matrix
regularMatrix <- acast(df, group.a ~ group.b, value.var = "estimated.difference.rate")
#max distance if NA
regularMatrix[is.na(regularMatrix)]<-1.0
distanceMatrix <- as.dist(regularMatrix)

#calculate dbscan
maxD<-max(distanceMatrix)
normD<-distanceMatrix/maxD
#no need to normalize - max_distance here is 1 and values we have are in that range
res <- dbscan(normD, eps = 1-similarity_threshold, minPts = 1)$cluster
names(res)<-labels(distanceMatrix)

# Format results
res.list <- lapply(split(res, names(res)), unname)
named_res <- lapply(res, function(x, prefix) paste0(prefix, x), prefix = "HaploGroup")
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
