#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--normal_files", default = NA, type = 'character', help = "normal samples input file names"),
				  make_option("--tumor_files", default = NA, type = 'character', help = "tumor samples input file names"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

in_file_normal = unlist(strsplit(x=opt$normal_files, split=" ", fixed=TRUE))
in_file_tumor = unlist(strsplit(x=opt$tumor_files, split=" ", fixed=TRUE))
out_file = opt$out_file

depth_n = list()
for (i in 1:length(in_file_normal)) {
	print(i)
	data = read.csv(file=in_file_normal[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% as.character(1:22)
	depth_n[[i]] = as.numeric(data[index,"depth"])
}
depth_n = do.call(cbind, depth_n)

depth_t = list()
for (i in 1:length(in_file_tumor)) {
	print(i)
	data = read.csv(file=in_file_tumor[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% as.character(1:22)
	depth_t[[i]] = as.numeric(data[index,"depth"])
}
depth_t = do.call(cbind, depth_t)

pca_n = prcomp(t(depth_n), center=TRUE, scale.=TRUE)
pca_t = predict(object=pca_n, newdata=t(depth_t))
x = c(pca_n$x[,1], pca_t[,1])
y = c(pca_n$x[,2], pca_t[,2])
bg = c(rep("grey90", nrow(pca_n$x)), rep("salmon", nrow(pca_t)))
col = c(rep("grey50", nrow(pca_n$x)), rep("black", nrow(pca_t)))
pch = 21

pdf(file=out_file, width=9, height=9)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(x=x, y=y, col = col, bg = bg, pch = pch, cex = 1, lwd = .1, axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "")
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "PC 1", line = 4, cex = 1.5)
mtext(side = 2, text = "PC 2", line = 4, cex = 1.5)
dev.off()
