#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(rjson)

get_region <- function(x) {

    y<-tail(unlist(strsplit(x, "/")),2)[1]

}

get_sample <- function(x) {

    y<-head(unlist(strsplit(x, "/")),3)[3]
    z<-head(unlist(strsplit(y, ".", fixed=TRUE)),1)

}

files<-sort(list.files(args[1], pattern="sorted_combos.tsv", full.names=T,recursive=T))
region<-unique(do.call(c,lapply(files, get_region)))
#all_list<-list()
tpr_list<-list()

for (r in region) {

    json_file<-file.path(args[1],"cluster", paste0(r, ".clusters.json"))
    clusters<-fromJSON(file=json_file)
    region_scores<-c()

    for (f in files[grepl(r,files)]) {
    
        df<-data.table::fread(f)
        sample<-get_sample(f)
        haps<-names(clusters)[grepl(sample,names(clusters))]
        true_clusters<-sort(c(clusters[[haps[1]]], clusters[[haps[2]]]))
        predicted_clusters<-sort(c(df$c1[1], df$c2[1]))
        region_scores<-c(region_scores, ifelse(identical(true_clusters,predicted_clusters),1,0))
        #all_list[[f]]<-data.frame(region=r, sample=sample, value=ifelse(identical(true_clusters,predicted_clusters),1,0))

    }

    tpr_list[[r]]<-rbind(c("TP",r,sum(region_scores)),c("FN",r,length(files[grepl(r,files)])-sum(region_scores)))


}

tpr_df<-do.call(rbind,tpr_list)
rownames(tpr_df)<-c()
df<-data.frame(ind=tpr_df[,1], variable=tpr_df[,2], value=as.numeric(tpr_df[,3]))
df$ind<-factor(df$ind, levels=c("FN","TP"))


#all_df<-do.call(rbind, all_list)
#rownames(all_df)<-c()

p<-ggplot(df, aes(x = variable, y = value, fill = ind)) + 
    geom_bar(position = "fill",stat = "identity", width=.5) +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=c("TP" = "darkred", "FN" = "darkblue"))

ggsave("tpr.pdf")