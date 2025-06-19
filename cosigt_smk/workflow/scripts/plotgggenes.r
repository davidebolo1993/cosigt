#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(gggenes)
library(data.table)
library(rjson)

setDTthreads(1)
gene_data <- fread(args[1])
cluster_mapping <- fromJSON(file=args[2])    
fai_index<-fread(args[3])

# add hapl start-end from fai
fai_index$start<-1
fai_index<-fai_index[,c(1,6,2)]
colnames(fai_index)<-c("molecule","mstart","mend")
gene_data<-merge(gene_data,fai_index,all.y=T)

# map clusters
clusters <- sapply(gene_data$molecule, function(mol) {
  cluster_mapping[[mol]]
})
gene_data$c <- clusters  

# simplify cluster names, so that they fit in the facets
gene_data$c <- gsub("HaploGroup", "", gene_data$c)

# order
gene_data <- gene_data[order(c, molecule)]
gene_data$start<-ifelse(is.na(gene_data$start),1,gene_data$start)
gene_data$end<-ifelse(is.na(gene_data$end),1,gene_data$end)
gene_data$strand<-ifelse(is.na(gene_data$strand),1,gene_data$strand)
#plot
p<-ggplot(gene_data, aes(xmin=start, xmax=end, y=molecule,
                         fill=gene, forward=strand, label=gene )) +
  geom_gene_arrow(aes(xmin=mstart, xmax=mend, y=molecule,forward = T),fill="#eeeeee",
                  color="#777777",show.legend = F,
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),arrow_body_height = unit(3,"mm"))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  show.legend=FALSE) +
  scale_fill_brewer(palette = "Set3") +
  geom_gene_label(align = "right") +
  theme_bw() +
  labs(y="haplotype") +
  facet_grid(rows=vars(c), scales="free_y", space="free_y") + 
  theme(
  ) +
  guides(fill = "none")+
  scale_y_discrete(
    expand = expansion(add = c(1,1)) # Add small expansion to both sides
  )

#save
plot_height <- max(5, 0.4*length(unique(gene_data$molecule)))
ggsave(args[4], height=plot_height, width=25, limitsize = FALSE)
