#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
load("modules/copy_number/CytoBand.RData")

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample", default = NA, type = 'character', help = "tumor sample"),
				 make_option("--type", default = NA, type = 'character', help = "type of plot"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'prunesegments.cn' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
			mcl = mean(x[indx,"Log2Ratio"])
			scl = sd(x[indx,"Log2Ratio"])
			ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
			x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
		} else {
			x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
		}
	}
	return(invisible(x))
}

if (opt$type=="raw") {

	infile = paste0("qdnaseq/bed/", opt$sample, ".bed")
	outfile = paste0("qdnaseq/copynumber/log2ratio/", opt$sample, ".pdf")
	data = read.table(file=infile, header=FALSE, sep="\t", skip=1, stringsAsFactors=FALSE)[,c(1,2,3,5),drop=FALSE]
	colnames(data) = c("Chromosome", "Start", "End", "Log2Ratio")
	col = rep("#9F6986", nrow(data))
	col[(data[,"Chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile, width=18, height=7)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"Log2Ratio"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2,2))
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
	axis(1, at = .5*(start+end), labels=labels, cex.axis = 1.15, las = 1)
	box(lwd=2.5)
	dev.off()
	
} else if (opt$type=="segmented") {

	infile = paste0("qdnaseq/copynumber/segmented/", opt$sample, ".RData")
	outfile = paste0("qdnaseq/copynumber/pcf/", opt$sample, ".pdf")
	load(infile)
	
	segmented = prunesegments.cn(x=segmented, n=7)
	lo2ratio = winsorize(data=data[,c("Chromosome","Start","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	end = NULL
	for (j in 1:22) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 22)
	start[2:22] = end[1:21]+1
	for (j in 1:22) {
		segmented[segmented[,"Chromosome"]==j,"Start"] = segmented[segmented[,"Chromosome"]==j,"Start"] + start[j]
		segmented[segmented[,"Chromosome"]==j,"End"] = segmented[segmented[,"Chromosome"]==j,"End"] + start[j]
		data[data[,"Chromosome"]==j,"Start"] = data[data[,"Chromosome"]==j,"Start"] + start[j]
	}
	col = "grey80"
	pdf(file=outfile, width=18, height=7)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"Start"], data[,"Log2Ratio"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2,2))
 	for (j in 1:nrow(segmented)) {
 		lines(x=c(segmented[j,"Start"], segmented[j,"End"]), y=rep(segmented[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
 	}
 	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:22) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(data[,"Start"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22), cex.axis = 0.85, las = 1)
	box(lwd=2.5)
	dev.off()

}
