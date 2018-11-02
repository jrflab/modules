#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

outfile_on_target = gsub(".cnr", ".on_target.pdf", opt$in_file, fixed=TRUE)
outfile_off_target = gsub(".cnr", ".off_target.pdf", opt$in_file, fixed=TRUE)

data = read.table(file=opt$in_file, header=TRUE, sep="\t", comment.char="#")

on_target = subset(data, data$gene=="-")
pdf(file=outfile_on_target, width=14)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(on_target$log2, type="n", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", ylim=c(-4,4))
points(on_target$log2, pch=19, cex=.5, col=ifelse(on_target$depth!=0, "steelblue", "red"))
abline(h=0, col="red", lty=1, lwd=1.5)
axis(1, at=NULL, cex.axis=1.5, padj=0.25)
axis(2, at=NULL, cex.axis=1.5, las=1)
mtext(side=1, text="Position (index)", line=4, cex=1.5)
mtext(side=2, text=expression(Log[2]~"Ratio"), line=4, cex=1.5)
box(lwd=2.5)
dev.off()

off_target = subset(data, data$gene!="-")
pdf(file=outfile_off_target, width=14)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(off_target$log2, type="n", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", ylim=c(-4,4))
points(off_target$log2, pch=19, cex=.5, col=ifelse(off_target$depth!=0, "steelblue", "red"))
axis(1, at=NULL, cex.axis=1.5, padj=0.25)
axis(2, at=NULL, cex.axis=1.5, las=1)
mtext(side=1, text="Position (index)", line=4, cex=1.5)
mtext(side=2, text=expression(Log[2]~"Ratio"), line=4, cex=1.5)
box(lwd=2.5)
dev.off()