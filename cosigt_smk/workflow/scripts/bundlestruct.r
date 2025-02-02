#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(tidyverse)
library(gggenes)
library(tidyverse)

# Set data.table threads
setDTthreads(1)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
bundles_in <- args[1]
bundles_struct <- args[2]
bundles_length <- args[3]
bundles_table <- args[4]

#read in and prepare for further operations
t <- fread(bundles_in, 
               header=FALSE, 
               sep="\t",
               col.names = c("contig",
                           "bstart",
                           "bend",
                           "binf"),
               skip = 1) %>%
                separate(binf, sep=":", 
                         c("bID",
                           "bdunno",
                            "bOrientation",
                            "entryix",
                            "exitix")) %>%
                mutate(bOrientation = as.numeric(bOrientation)) %>%
                group_by(contig) %>%
                mutate(seg_len = bend-bstart) %>%
                mutate(rel_end = cumsum(seg_len)) %>%
                mutate(rel_start = lag(rel_end,default=0)) %>%
                separate(contig,c("sample_id","sample_hap", "sample_contig"),sep="#|\\.",remove=FALSE) %>%
                dplyr::select(-c(sample_hap,sample_contig)) #remove those 2 - don't needed

#bundles struct
hap_sums <- t %>% 
            mutate(bID_d = paste(bID,".",bOrientation,sep="")) %>%
            group_by(contig) %>%
            summarize(hap_struc = paste0(bID,collapse="-"),
                      hap_struc_d = paste0(bID_d,collapse="-"))

contig_structs <- hap_sums %>% 
  dplyr::select(contig, hap_struc_d)

#write to file
fwrite(contig_structs,bundles_struct,sep="\t",row.names = FALSE,quote=FALSE)

#now prepare to calculate lengths
t <- inner_join(t,hap_sums,by="contig")

hap_info <- t %>% group_by(contig,sample_id,hap_struc,hap_struc_d) %>%
                 summarize(hap_block_start = min(bstart))

t<-inner_join(t,hap_info,by=c("contig","hap_struc", "hap_struc_d", "sample_id")) %>% 
  mutate(rel_start2 = bstart-hap_block_start,
  rel_end2 = bend-hap_block_start)

fwrite(t, bundles_table, sep="\t",row.names = FALSE,quote=FALSE)

#get representative structs
representative_haps <- t %>% group_by(contig) %>%
  filter(row_number()==1) %>% 
  ungroup() %>%
  group_by(hap_struc_d) %>%
  mutate(n = n()) %>%
  filter(row_number()==1) %>% 
  ungroup() %>%
  dplyr::select(contig,n,hap_struc_d)

#get bundle lengths
bID_lens <- t %>% filter(contig %in% representative_haps$contig) %>% 
  mutate(bID_d = paste(bID,bOrientation,sep=".")) %>%
  mutate(l = bend-bstart) %>%
  group_by(bID_d) %>%
  summarize(l = median(l))

#write to file
fwrite(bID_lens,bundles_length,sep="\t",row.names = FALSE,quote=FALSE)

#plot structures
g<-ggplot(t) +
  geom_gene_arrow(aes(xmin=rel_start2,
    xmax=rel_end2,
    y=contig,
    fill=bID,
    forward=!bOrientation),
    arrowhead_height = unit(1, "mm"), 
    arrowhead_width = unit(1, "mm"),
    arrow_body_height = unit(.5, "mm"),
    position=position_nudge(x=0,y=-.2),
    size=.1)+
  geom_text(aes(x=rel_start2,
    y=contig,
    label=bID),
    size=0) +
  theme_bw(base_size=4)+
  theme(legend.position="None")

#save graphical output
ggsave(gsub("tsv", "pdf", bundles_struct), width=15, height=20)
