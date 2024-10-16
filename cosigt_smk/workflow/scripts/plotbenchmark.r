#!/usr/bin/Rscript
library(data.table)
library(ggplot2)

#do not exceed
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)

get_method <- function(x) {
    y<-tail(unlist(strsplit(x, ".", fixed=T)),3)[1]
}

get_region <- function(x) {
    y<-tail(unlist(strsplit(x, ".", fixed=T)),4)[1]
    z<-tail(unlist(strsplit(y, "/")),1)
}

#input files, divided by type
files<-sort(list.files(args[1], pattern=".benchmark.txt", full.names=T))
#those having a region in their name
files_wregion<-grep("^[^_]+_[^_]+_[^_]+_[^_]+[^.]_", files, value = TRUE)
#those that do not have a region in their name
files_noregion<-grep("^[^_]+_[^_]+_[^_]+_[^_]+[^.]_", files, invert=TRUE, value = TRUE)

#iterate over methods and aggregate where we have many
method<-unique(do.call(c,lapply(files_wregion, get_method)))
region<-unique(do.call(c,lapply(files_wregion, get_region)))
benchmark_list<-list()

l<-1
for (m in method) {
    for (r in region) {
        f1<-files_wregion[grep(m,files_wregion)]
        f2<-f1[grep(r,f1)]
        values_time<-c()
        values_mem<-c()
        for (i in c(1:length(f2))) {
            df<-fread(f2[i], header=T)
            values_time<-c(values_time,df$s)
            values_mem<-c(values_mem,df$max_rss)
        }
        benchmark_list[[l]]<-data.frame(region=r,method=m,value.time=values_time, value.mem=values_mem)
        l<-l+1
    }   
}

#add rules that are run just once in total
for (f in files_noregion) {
    l<-l+1
    x<-tail(unlist(strsplit(f, ".", fixed=T)),4)[1]
    method<-tail(unlist(strsplit(x, "/")),1)
    df<-fread(f,header=T)
    benchmark_list[[l]]<-data.frame(region="None",method=method,value.time=df$s, value.mem=df$max_rss)
}

benchmark_df<-do.call(rbind,benchmark_list)

p<-ggplot(benchmark_df, aes(x=method,y=value.time, fill=region)) + 
    geom_boxplot() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="top")+
    scale_y_log10(breaks=c(.1,1,10,100,1000)) +
    labs(x="step", y="runtime (s)") 

ggsave(args[2], width=20)

p<-ggplot(benchmark_df, aes(x=method,y=value.mem, fill=region)) + 
    geom_boxplot() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="top")+
    scale_y_log10(breaks=c(.1,1,10,100,1000,10000,100000)) +
    labs(x="step", y="max_rss (MB)")

ggsave(args[3], width=20)