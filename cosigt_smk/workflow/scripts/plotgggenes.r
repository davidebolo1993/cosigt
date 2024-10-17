#!/usr/bin/Rscript
## TO BE STRUCTURED AGAIN

args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(gggenes)
library(rjson)

x <- data.table::fread(args[1])
y <- fromJSON(file=args[2])
reference<-args[3]

if (!(reference %in% x$molecule)) {

    x$strand<-abs(x$strand-1)

}

x$old_label<-x$molecule
x$molecule<-gsub("_inv", "", x$molecule)
c<-c()

for (i in c(1:nrow(x))) {

    if (x$molecule[i] %in% names(y)) {

        c<-c(c, y[[x$molecule[i]]])
    
    } else {

        c<-c(c, "reference")
    }

}

x$c<-c
x<-x[order(c),]
x<-subset(x, c != "reference")

x$ymin<-0
x$ymax<-0
x$molecule<-as.numeric(factor(x$molecule, levels=unique(x$molecule)))

for (cl in unique(x$c)) {

    s<-subset(x, (c == cl))
    first<-as.numeric(head(s,1)$molecule)-0.4
    last<-as.numeric(tail(s,1)$molecule)+0.4
    x[which(x$c == cl),]$ymin<-first
    x[which(x$c == cl),]$ymax<-last

}

p<-ggplot(x, aes(xmin=start, xmax=end, y=molecule, fill=gene, forward=strand, label=gene)) + 
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), show.legend=FALSE) + 
    scale_fill_brewer(palette = "Set3") +
    geom_rect(data= data.frame(ymin=unique(x$ymin),ymax=unique(x$ymax),cluster=unique(x$c)), inherit.aes = FALSE,  mapping=aes(xmin = -Inf, xmax = +Inf, ymin = ymin, ymax = ymax, color = cluster), alpha=0.4, show.legend=TRUE, linewidth=1, fill=NA) +
    #scale_color_brewer(palette = "Dark2") +
    geom_gene_label(align = "right") +
    theme_genes() +
    labs(y="haplotype") +
    scale_y_continuous(breaks=c(1:length(unique(x$old_label))), labels=unique(x$old_label))+
    guides(fill = "none")

ggsave(paste0(args[3], ".pdf"), height=20, width=20)