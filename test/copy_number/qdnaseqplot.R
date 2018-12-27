#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample", default = NA, type = 'character', help = "tumor sample"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

infile = paste0("qdnaseq/bed/", opt$sample, ".bed")
outfile = paste0("qdnaseq/copynumber/log2ratio/", opt$sample, ".pdf")

data = read.table(file=infile, header=FALSE, sep="\t", skip=1, stringsAsFactors=FALSE)[,c(1,2,3,5),drop=FALSE]
colnames(data) = c("Chromosome", "Start", "End", "Log2Ratio")
col = rep("#9F6986", nrow(data))
col[(data[,"Chromosome"]%%2)==1] = "#CECAC5"
pdf(file=outfile, width=18, height=7)
par(mar=c(5, 5, 4, 2)+.1)
plot(data[,"Log2"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
abline(v=1, col="goldenrod3")
abline(h=0, col="red")
for (j in 2:max(data[,"Chromosome"])) {
v = min(which(data[,"Chromosome"]==j))
	abline(v=v, col="goldenrod3")
}
abline(v=max(nrow(data)), col="goldenrod3")
start = NULL
end = NULL
for (j in 1:max(data[,"Chromosome"])) {
	start[j] = min(which(data[,"Chromosome"]==j))
	end[j] = max(which(data[,"Chromosome"]==j))
}
labels = 1:max(data[,"Chromosome"])
labels[labels==23] = "X"
axis(1, at = .5*(start+end), labels=labels, cex.axis = 0.85, las = 1)
box(lwd=2.5)
dev.off()
