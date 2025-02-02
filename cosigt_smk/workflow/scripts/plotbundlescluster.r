#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(gggenes)
library(data.table)
library(rjson)

# Set data.table threads
setDTthreads(1)

bundles_data <- fread(args[1])
cluster_mapping <- fromJSON(file=args[2])
        
# Map clusters
clusters <- sapply(bundles_data$contig, function(contig) {
    return(cluster_mapping[[contig]])
})
    
bundles_data$c <- clusters
    
# Order data
bundles_data <- bundles_data[order(clusters)]

#calculate data with cluster boundary information
bundles_data$ymin <- 0
bundles_data$ymax <- 0
bundles_data$old_label <- bundles_data$contig
bundles_data$contig <- as.numeric(factor(bundles_data$contig, levels=unique(bundles_data$contig)))
    
# Calculate boundaries for each cluster
for (cl in unique(bundles_data$c)) {
    subset_data <- subset(bundles_data, (c == cl))
    first <- as.numeric(head(subset_data, 1)$contig) - 0.4
    last <- as.numeric(tail(subset_data, 1)$contig) + 0.4
    bundles_data[which(bundles_data$c == cl),]$ymin <- first
    bundles_data[which(bundles_data$c == cl),]$ymax <- last
}

#plot
g<-ggplot(bundles_data,aes(xmin=rel_start2,
    xmax=rel_end2,
    y=contig,
    fill=bID,
    forward=!bOrientation, label=bID)) +
  geom_gene_arrow(,
    arrowhead_height = unit(1, "mm"), 
    arrowhead_width = unit(1, "mm"),
    arrow_body_height = unit(.5, "mm"),
    position=position_nudge(x=0,y=-.2),
    size=.1)+
    geom_rect(data = data.frame(
        ymin=unique(bundles_data$ymin),
        ymax=unique(bundles_data$ymax),
        cluster=unique(bundles_data$c)
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
        breaks=seq_len(length(unique(bundles_data$old_label))), 
        labels=unique(bundles_data$old_label)
    ) +
    guides(fill = "none")

#save
ggsave(args[3], height=20, width=20)
