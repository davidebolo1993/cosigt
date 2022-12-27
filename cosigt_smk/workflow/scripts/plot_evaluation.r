args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)

#binary evaluation
binev<-function(df,items) {
  binary_evaluation<-rep(0,nrow(df))
  for (i in c(1:nrow(df))) {
    if (df$V2[i] == 2) {
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
outpdf<-file.path(args[3])

#list.files
files<-list.files(infolder, recursive=T, pattern="evaluation.tsv", full.names=T)

all_files<-list()
c<-0

for (f in files) {
  
  c<-c+1
  tab<-read.table(f, header=F, sep="\t")
  tab<-head(tab,items)
  tab$V4<-binev(tab,items)
  all_files[[c]]<-tab

}

#tables converted to a binary evaluation:
#1. indicates that the true combination has been identified
#0. indicates that the true combination has not been identified yet
df<-data.table(do.call(rbind,all_files))

result<-rep(0,items)
#sumarise across samples
for (i in c(1:items)) {
  sub_df<-df[V1==i]
  result[i] = sum(sub_df$V4)/length(unique(sub_df$V3))
}

#plot
plot_df<-data.frame(x=c(1:items), y=result)
p<-ggplot(plot_df, aes(x=x,y=y)) + geom_point() + geom_line() +
  scale_y_continuous(limits = c(0.0,1.0)) + 
  scale_x_continuous(limits=c(1,items), n.breaks = items) +
  theme_bw() +
  xlab("#rank") + ylab("TPR")

#save
ggsave(p, filename = outpdf)
