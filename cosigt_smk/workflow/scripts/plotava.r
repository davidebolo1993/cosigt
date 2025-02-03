#!/usr/bin/Rscript
#IS THIS SOMEHOW USEFUL?

args <- commandArgs(trailingOnly = TRUE)
library(SVbyEye)

input_paf<-args[1]
output_pdf<-args[2]
ref_path<-args[3]

#read
paf.table <- readPaf(paf.file = input_paf, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#subset to single seq vs target and plot
sub.paf<-subset(paf.table, (grepl(ref_path, t.name) & !grepl(ref_path,q.name)))
seqnames<-unique(sub.paf$q.name)
pltlist<-list()
for (i in c(1:length(seqnames))) {
        simplify<-paste(unlist(strsplit(seqnames[i], "#"))[c(1,2)],collapse="#")
        sub.sub.paf<-subset(sub.paf, (q.name == seqnames[i]))
        pltlist[[simplify]]<-plotMiro(paf.table = sub.sub.paf, binsize = 1000)
}
for (l in c(1:length(pltlist))) {
        pdf(file.path(dirname(output_pdf),paste0(names(pltlist)[l], "_to_", ref_path, ".pdf")), width=20, height=5)
        pltlist[[l]]
        dev.off()
}
#plot all vs all
pdf(output_pdf, width=20, height=10)
plotAVA(paf.table = paf.table, binsize=5000)
dev.off()