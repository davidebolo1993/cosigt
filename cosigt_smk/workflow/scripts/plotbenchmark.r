#!/usr/bin/Rscript

library(ggplot2)


get_method <- function(x) {

    y<-tail(unlist(strsplit(x, ".", fixed=T)),3)[1]

}

get_region <- function(x) {

    y<-tail(unlist(strsplit(x, ".", fixed=T)),4)[1]
    z<-tail(unlist(strsplit(y, "/")),1)

}

files<-sort(list.files("benchmarks", pattern=".benchmark.txt", full.names=T))
method<-unique(do.call(c,lapply(files, get_method)))
region<-unique(do.call(c,lapply(files, get_region)))
benchmark_list<-list()

l<-1

for (m in method) {
    
    for (r in region) {

        f1<-files[grep(m,files)]
        f2<-f1[grep(r,f1)]
        values<-c()

        for (i in c(1:length(f2))) {

            df<-data.table::fread(f2[i], header=T)
            values<-c(values,df$s)

        }

        benchmark_list[[l]]<-data.frame(region=r,method=m,values=values)
        l<-l+1

    }   

}


benchmark_df<-do.call(rbind,benchmark_list)
p<-ggplot(benchmark_df, aes(x=method,y=values, fill=region)) + 
    geom_boxplot() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="top")+
    scale_y_log10(breaks=c(.1,1,10,100,1000)) +
    labs(x="step", y="runtime (s)")
    

ggsave("benchmark.pdf", p, width=20)
