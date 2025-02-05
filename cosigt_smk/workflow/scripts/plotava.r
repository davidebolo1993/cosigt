#!/usr/bin/Rscript
#STILL DON'T KNOW IF THIS IS FULLY USEFUL

args <- commandArgs(trailingOnly = TRUE)
library(SVbyEye)

input_paf<-args[1]
output_pdf<-args[2]
ref_path<-args[3]

#read
paf.table <- readPaf(paf.file = input_paf, include.paf.tags = TRUE, restrict.paf.tags = "cg")
paf.table<-flipPaf(paf.table, majority.strand="+")
#subset to single seq vs target and plot
sub.paf<-subset(paf.table, (grepl(ref_path, t.name) & !grepl(ref_path,q.name)))
seqnames<-unique(sub.paf$q.name)
for (i in 1:length(seqnames)) {
    simplify<-paste(unlist(strsplit(seqnames[i], "#"))[c(1,2)],collapse="#")
    sub.sub.paf<-subset(sub.paf, (q.name == seqnames[i]))
    pdf(file.path(dirname(output_pdf),paste0(simplify, "_to_", ref_path, ".pdf")), width=20, height=5)
    print(plotMiro(paf.table = sub.sub.paf, binsize = 1000))
    dev.off()
}
#plot all vs all
pdf(output_pdf, width=20, height=10)
plotAVA(paf.table = paf.table, binsize=5000)
dev.off()
