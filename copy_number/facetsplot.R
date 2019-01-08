#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load(opt$in_file)
data = out2$jointseg
data[data[,"chrom"]=="X", "chrom"] = 23
data[data[,"chrom"]=="Y", "chrom"] = 24
data[,"chrom"] = as.numeric(data[,"chrom"])
data = subset(data, data[,"chrom"]<=23)
tmp = data[,c("chrom", "maploc", "cnlr"),drop=FALSE]
colnames(tmp) = c("Chromosome", "Position", "Log2Ratio")
tmp = winsorize(data=tmp, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
data[tmp$wins.outliers[,3]!=0,"cnlr"] = NA
col = rep("#9F6986", nrow(data))
col[(data[,"chrom"]%%2)==1] = "#CECAC5"
pdf(file=opt$out_file, width=14, height=5)
par(mar=c(5, 5, 4, 2)+.1)
plot(data[,"cnlr"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
abline(v=1, col="goldenrod3")
abline(h=0, col="red")
for (j in 2:max(data[,"chrom"])) {
	v = min(which(data[,"chrom"]==j))
	abline(v=v, col="goldenrod3")
}
abline(v=max(nrow(data)), col="goldenrod3")
start = NULL
end = NULL
for (j in 1:max(data[,"chrom"])) {
	start[j] = min(which(data[,"chrom"]==j))
	end[j] = max(which(data[,"chrom"]==j))
}
labels = 1:max(data[,"chrom"])
labels[labels==23] = "X"
axis(1, at = .5*(start+end), labels=labels, cex.axis = 0.85, las = 1)
box(lwd=2.5)
dev.off()

