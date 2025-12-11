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
# added by chiara.paleni (if the .medoids.tsv is in args[5])
medoids<-fread(args[5],header=F)

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

# order
gene_data <- gene_data[order(c, molecule)]
gene_data$start<-ifelse(is.na(gene_data$start),1,gene_data$start)
gene_data$end<-ifelse(is.na(gene_data$end),1,gene_data$end)
gene_data$strand<-ifelse(is.na(gene_data$strand),1,gene_data$strand)
#plot
p <- ggplot(gene_data, aes(xmin=start, xmax=end, y=molecule,
                           fill=gene, forward=strand, label=gene)) +
  geom_gene_arrow(aes(xmin=mstart, xmax=mend, y=molecule, forward = TRUE),
                  fill="#eeeeee", color="#777777", show.legend = FALSE,
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),
                  arrow_body_height = unit(3, "mm")) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  show.legend = TRUE) +  # allow legend here
  scale_fill_brewer(palette = "Set3", name = "Gene") +  # legend title
  geom_gene_label(align = "right") +
  theme_bw() +
  labs(y = "Haplotypes", x = "Genomic Position") +
  facet_grid(rows=vars(c), scales = "free_y", space = "free_y") +
  theme(strip.text.y.right = element_text(angle = 0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.6, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)
  ) +
  scale_y_discrete(expand = expansion(add = c(1, 1)))

#save
plot_height <- max(5, 0.4*length(unique(gene_data$molecule)))
ggsave(args[4], height=plot_height, width=25, limitsize = FALSE)

#chiara.paleni added this
mini_gene<-gene_data[gene_data$molecule %in% medoids$V2,]
unique_mol<-gene_data[!duplicated(gene_data$molecule),]
abundance<-data.frame(table(unique_mol$c))
mini_gene<-merge(mini_gene,abundance,by.x="c",by.y="Var1")
mini_gene$label<-paste0(mini_gene$c,"\n",mini_gene$Freq," haplotype(s)")
mini_gene<-mini_gene[order(mini_gene$Freq,decreasing = T),]
mini_gene$label<-factor(mini_gene$label,levels=unique(mini_gene$label))

p <- ggplot(mini_gene, aes(xmin=start, xmax=end, y=molecule,
                           fill=gene, forward=strand, label=gene)) +
  geom_gene_arrow(aes(xmin=mstart, xmax=mend, y=molecule, forward = TRUE),
                  fill="#eeeeee", color="#777777", show.legend = FALSE,
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),
                  arrow_body_height = unit(3, "mm")) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  show.legend = TRUE) +  # allow legend here
  scale_fill_brewer(palette = "Set3", name = "Gene") +  # legend title
  geom_gene_label(align = "right") +
  theme_bw() +
  labs(y = "Haplotypes", x = "Genomic Position") +
  facet_grid(rows=vars(label), scales = "free_y", space = "free_y") +
  theme(strip.text.y.right = element_text(angle = 0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.6, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)
  ) +
  scale_y_discrete(expand = expansion(add = c(1, 1)))

plot_height <- max(5, 0.4*length(unique(medoids$V1)))
ggsave(gsub(".png",".representative.png",args[4]), height=plot_height, width=25, limitsize = FALSE)