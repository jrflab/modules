#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

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

bin_size = as.numeric(data[index,"end"]) - as.numeric(data[index,"start"])
var_bin_n = apply(depth_n, 1, sd, na.rm=TRUE)
var_bin_t = apply(depth_t, 1, sd, na.rm=TRUE)
data = data.frame(bin_size, var_bin_n, var_bin_t)
write.table(data, file=out_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

ymin = min(var_bin_n, var_bin_t)
ymax = max(var_bin_n, var_bin_t)

pdf(file=gsub(".tsv", ".pdf", x=out_file, fixed=TRUE), width=7, height=7)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(bin_size, var_bin_n, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", log="y", ylim=c(ymin, ymax))
points(x=bin_size, y=var_bin_n, col = "grey50", bg = "grey90", pch = 21, cex = 1, lwd = .1)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "Bin size (bp)", line = 4, cex = 1.5)
mtext(side = 2, text = "SD", line = 5, cex = 1.5)
plot(bin_size, var_bin_t, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", log="y", ylim=c(ymin, ymax))
points(x=bin_size, y=var_bin_t, col = "black", bg = "steelblue", pch = 21, cex = 1, lwd = .1)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "Bin size (bp)", line = 4, cex = 1.5)
mtext(side = 2, text = "SD", line = 5, cex = 1.5)
plot(var_bin_n, var_bin_t, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", log="xy", xlim=c(ymin, ymax), ylim=c(ymin, ymax))
points(x=var_bin_n, y=var_bin_t, col = "black", bg = "steelblue", pch = 21, cex = 1, lwd = .1)
abline(a=0, b=1, col="goldenrod3", lwd=2)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "Normal SD", line = 4, cex = 1.5)
mtext(side = 2, text = "Tumor SD", line = 5, cex = 1.5)
dev.off()
