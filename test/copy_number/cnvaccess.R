#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {
	
	dataA = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".A.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataA[dataA[,"chromosome"]=="X", "chromosome"] = 23
	dataA[dataA[,"chromosome"]=="Y", "chromosome"] = 24
	dataA[,"chromosome"] = as.numeric(dataA[,"chromosome"])
	dataA = subset(dataA, dataA[,"chromosome"]<=23)
	
	dataB = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".B.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataB[dataB[,"chromosome"]=="X", "chromosome"] = 23
	dataB[dataB[,"chromosome"]=="Y", "chromosome"] = 24
	dataB[,"chromosome"] = as.numeric(dataB[,"chromosome"])
	dataB = subset(dataB, dataB[,"chromosome"]<=23)
	
	dataC = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".C.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataC[dataC[,"chromosome"]=="X", "chromosome"] = 23
	dataC[dataC[,"chromosome"]=="Y", "chromosome"] = 24
	dataC[,"chromosome"] = as.numeric(dataC[,"chromosome"])
	dataC = subset(dataC, dataC[,"chromosome"]<=23)
	
	data = rbind(dataA, dataB, dataC)
	index = order(data[,2])
	data = data[index,,drop=FALSE]
	index = order(data[,1])
	data = data[index,,drop=FALSE]
	
	pdf(file=outfile_f = paste0("cnvaccess/plot/log2/", opt$sample_name, ".pdf"), width=10, height=4.2)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"log2"], type="p", pch=19, cex=.25, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:max(data[,"chromosome"])) {
		v = min(which(data[,"chromosome"]==j))
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(nrow(data)), col="goldenrod3")
	start = NULL
	end = NULL
	for (j in 1:max(data[,"chromosome"])) {
		start[j] = min(which(data[,"chromosome"]==j))
		end[j] = max(which(data[,"chromosome"]==j))
	}
	labels = 1:max(data[,"chromosome"])
	labels[labels==23] = "X"
	axis(1, at = .5*(start+end), labels=labels, cex.axis = 0.85, las = 1)
	box(lwd=1.5)
	dev.off()
	
}
