#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(gggenes)
library(data.table)
library(rjson)

# Set data.table threads
setDTthreads(1)

gene_data <- fread(args[1])
cluster_mapping <- fromJSON(file=args[2])
    
# Clean molecule names
gene_data$molecule <- gsub("_inv", "", gene_data$molecule)
    
# Map clusters
clusters <- sapply(gene_data$molecule, function(mol) {
    if (mol %in% names(cluster_mapping)) {
        return(cluster_mapping[[mol]])
    }
    return("reference")
})
    
gene_data$c <- clusters
    
# Filter and order data
gene_data <- gene_data[order(clusters)]
gene_data <- subset(gene_data, clusters != "reference")


#calculate data with cluster boundary information
gene_data$ymin <- 0
gene_data$ymax <- 0
gene_data$old_label <- gene_data$molecule
gene_data$molecule <- as.numeric(factor(gene_data$molecule, levels=unique(gene_data$molecule)))
    
# Calculate boundaries for each cluster
for (cl in unique(gene_data$c)) {
    subset_data <- subset(gene_data, (c == cl))
    first <- as.numeric(head(subset_data, 1)$molecule) - 0.4
    last <- as.numeric(tail(subset_data, 1)$molecule) + 0.4
    gene_data[which(gene_data$c == cl),]$ymin <- first
    gene_data[which(gene_data$c == cl),]$ymax <- last
}


#plot
p<-ggplot(gene_data, aes(xmin=start, xmax=end, y=molecule, 
                        fill=gene, forward=strand, label=gene)) + 
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
        arrowhead_width = unit(1, "mm"), 
        show.legend=FALSE) + 
    scale_fill_brewer(palette = "Set3") +
    geom_rect(data = data.frame(
        ymin=unique(gene_data$ymin),
        ymax=unique(gene_data$ymax),
        cluster=unique(gene_data$c)
    ), 
    inherit.aes = FALSE,  
    mapping=aes(xmin = -Inf, xmax = +Inf, 
                ymin = ymin, ymax = ymax, 
                color = cluster), 
        alpha=0.4, show.legend=TRUE, 
        linewidth=1, fill=NA) +
    geom_gene_label(align = "right") +
    theme_genes() +
    labs(y="haplotype") +
    scale_y_continuous(
        breaks=seq_len(length(unique(gene_data$old_label))), 
        labels=unique(gene_data$old_label)
    ) +
    guides(fill = "none")


ggsave(args[3], height=20, width=20)