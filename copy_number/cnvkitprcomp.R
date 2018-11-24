#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file names"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

col = "steelblue"
pch = 1
cex = 1

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
pca = prcomp(t(depth), center=TRUE, scale.=TRUE)
pdf(file=out_file, width=9, height=9)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(x=pca$x[,1], y=pca$x[,2], col = col, pch = pch, cex = cex, axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "")
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25)
axis(2, at = NULL, cex.axis = 1.5, las = 1)
mtext(side = 1, text = "PC 1", line = 4, cex = 1.5)
mtext(side = 2, text = "PC 2", line = 4, cex = 1.5)
dev.off()
