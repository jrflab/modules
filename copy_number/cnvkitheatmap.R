#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file names"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

in_file = unlist(strsplit(x=opt$in_file, split=" ", fixed=TRUE))
out_file = opt$out_file

depth = list()
for (i in 1:length(in_file)) {
	print(i)
	data = read.csv(file=in_file[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% c(as.character(1:22), "X")
	depth[[i]] = as.numeric(data[index,"depth"])
}
depth = do.call(cbind, depth)
pdf(file=out_file, width=14, height=14)
heatmap(x=depth, labRow=rep(" ", nrow(depth)), labCol=rep(" ", ncol(depth)), col=colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
dev.off()

png(file=gsub(".pdf", ".png", out_file, fixed=TRUE), width=1440, height=1440)
heatmap(x=depth, labRow=rep(" ", nrow(depth)), labCol=rep(" ", ncol(depth)), col=colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
dev.off()
