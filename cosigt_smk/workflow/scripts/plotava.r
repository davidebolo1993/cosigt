#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
library(SVbyEye)

# move saving to function - this should help understanind why sometimes we miss a plot even if no error in smk?
save_plot_as_png <- function(plot_obj, file_path, width = 20, height = 5, res = 300) {
  tryCatch({
    png(file_path, width = width, height = height, units = "in", res = res)
    print(plot_obj)
  }, error = function(e) {
    message("An error occurred while creating the plot: ", e$message)
  }, finally = {
    # ensure the device is closed
    if (dev.cur() != 1) {
      dev.off()
    }
  })
}

# args in
input_paf <- args[1]
output_png <- args[2]
ref_path <- args[3]

# read paf
paf.table <- readPaf(paf.file = input_paf, include.paf.tags = TRUE, restrict.paf.tags = "cg")

# single sequence vs target
sub.paf <- subset(paf.table, (grepl(ref_path, t.name) & !grepl(ref_path, q.name)))
seqnames <- unique(sub.paf$q.name)

if (length(seqnames) > 0) {
  
  #individual plots
  for (seqname in seqnames) {
    simplify_name <- paste(unlist(strsplit(seqname, "#"))[c(1, 2)], collapse = "#")
    sub.sub.paf <- subset(sub.paf, q.name == seqname)

    if (sum(sub.sub.paf$strand == "-") > sum(sub.sub.paf$strand == "+")) {
      sub.sub.paf <- flipPaf(paf.table = sub.sub.paf, flip.seqnames = seqname)
      paf.table <- flipPaf(paf.table = paf.table, flip.seqnames = seqname)
    }

    miro_plot <- plotMiro(paf.table = sub.sub.paf, binsize = 1000)
    output_file <- file.path(dirname(output_png), paste0(simplify_name, "_to_", ref_path, ".png"))
    save_plot_as_png(miro_plot, output_file)
  }

  # all-vs-all
  ava_plot <- plotAVA(paf.table = paf.table, binsize = 1000)
  save_plot_as_png(ava_plot, output_png, height = 10)

} else {
  
  #create empty in case
  file.create(output_png)
  message("No sequences found to plot. Created an empty file.")
}
