args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(ggforce)
library(rjson)

#binary evaluation
binev<-function(df,maxtimes,cluster) {

  binary_evaluation<-rep(0,nrow(df))
  cluster_to_consider<-grep(df$V3[1],unique(c(df$V4,df$V5)),value=TRUE) #find elements we need to consider

  all_1<-c(cluster_to_consider[1])

  for (l in cluster[[cluster_to_consider[1]]]) {

    all_1<-append(all_1, l)

  }

  if (maxtimes == 2) {

    all_2<-c(cluster_to_consider[2])

    for (l in cluster[[cluster_to_consider[2]]]) {

      all_2<-append(all_2, l)

    }

  } else { #just one haplotype in the graph, second can be anything

    all_2<-grep(df$V3[1],unique(c(df$V4,df$V5)),value=TRUE,invert=TRUE) #can be anything

  }

  for (i in c(1:nrow(df))) {
        
    if ((df$V4[i] %in% all_1 && df$V5[i] %in% all_2) || (df$V4[i] %in% all_2 && df$V5[i] %in% all_1)) { #must be this way here

      binary_evaluation[c(i:nrow(df))]<-rep(1,nrow(df)-i+1)
      break

    } else {

      binary_evaluation[i] == 0

    }
  }
  
  return(binary_evaluation)

}

#inputs
infolder<-file.path(args[1])
items<-as.integer(args[2])

if (args[3] != "") {

  cluster<-fromJSON(file = file.path(args[3]))

}

outpdf<-file.path(args[4])

#list.files with evaluation
files<-list.files(infolder, recursive=T, pattern="evaluation.tsv", full.names=T)

#list.files with max number of times
maxtimes<-list.files(infolder, recursive=T, pattern="maxtimes.txt", full.names=T)

all_files<-list()

for (i in c(1:length(files))) {
  
  tab<-fread(files[i], header=F, sep="\t")
  cmp<-fread(maxtimes[i], header=F,sep="\t")

  if (args[3] == "") {

    cluster<-list()
    id<-basename(dirname(files[i]))
    vals<-unique(grep(id, c(tab$V4,tab$V5), value=TRUE))

    for (v in vals) {

      cluster[[v]]<-list()

    }

  }

  tab$V7<-binev(tab, as.integer(cmp$V1[1]), cluster)
  tab$V8 <- ifelse(tab$V1 <= items,TRUE,FALSE)
  all_files[[i]]<-tab

}

#tables converted to a binary evaluation:
#1. indicates that the true combination has been identified
#0. indicates that the true combination has not been identified yet
df<-data.table(do.call(rbind,all_files))

result<-rep(0,nrow(tab))
#sumarise across samples
for (i in c(1:length(result))) {
  sub_df<-df[V1==i]
  result[i] = sum(sub_df$V7)/length(unique(sub_df$V3))
}

#plot
plot_df<-data.frame(x=c(1:length(result)), y=result)
p<-ggplot(plot_df, aes(x=x,y=y)) + geom_point() + geom_line() +
  scale_y_continuous(limits = c(0.0,1.0)) + 
  scale_x_continuous(limits=c(1,length(result))) +
  theme_bw() +
  xlab("#rank") + ylab("TPR") +
  facet_zoom(xlim=c(0,items),ylim=c(0,1), horizontal=FALSE, zoom.data=V8)

#save
ggsave(p, filename = outpdf)
