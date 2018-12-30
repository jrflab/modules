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
				 make_option("--type", default = NA, type = 'character', help = "type of plot"),
				 make_option("--rho", default = NA, type = 'numeric', help = "tumor purity"),
				 make_option("--psi", default = NA, type = 'numeric', help = "tumor ploidy"),
				 make_option("--gamma", default = NA, type = 'numeric', help = "log2 ratio compression"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

opt$rho = ifelse(is.na(as.numeric(opt$rho)), 1, as.numeric(opt$rho))
opt$psi = ifelse(is.na(as.numeric(opt$psi)), 2, as.numeric(opt$psi))
opt$gamma = ifelse(is.na(as.numeric(opt$gamma)), 1, as.numeric(opt$gamma))

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
	axis(1, at = .5*(start+end), labels=c(1:22), cex.axis = 1.15, las = 1)
	for (k in 1:8) {
		abline(h=(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), col="darkorange", lty=3)
		mtext(text=k, side=4, line=.5, at=(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), las=2, cex=.75, col="darkorange")
	}
	box(lwd=2.5)
	dev.off()

} else if (opt$type=="bychromosome") {

	infile = paste0("qdnaseq/copynumber/segmented/", opt$sample, ".RData")
	if (!dir.exists("qdnaseq/copynumber/bychr/")) {
		dir.create("qdnaseq/copynumber/bychr/")
	}
	if (!dir.exists(paste0("qdnaseq/copynumber/bychr/", opt$sample, "/"))) {
		dir.create(paste0("qdnaseq/copynumber/bychr/", opt$sample, "/"))
	}
	load(infile)
	segmented = prunesegments.cn(x=segmented, n=7)
	for (ii in 1:22) {
		pdf(file=paste0("qdnaseq/copynumber/bychr/", opt$sample, "/", ii, ".pdf"))
		zz = split.screen(figs=matrix(c(0,1,.15,1, 0.065,.975,0.1,.4), nrow=2, ncol=4, byrow=TRUE))
		screen(zz[1])
		par(mar = c(6.1, 6, 4.1, 3))
		start = 1
		end = max(CytoBand[CytoBand[,"Chromosome"]==ii,"End"])
		plot(1, 1, type="n", xlim=c(start,end), ylim=c(-2,2), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
		index = data[,"Chromosome"]==ii
		points(data[index,"Start"], data[index,"Log2Ratio"], type="p", pch=".", cex=1.15, col="grey80")
		tmp = subset(segmented, segmented[,"Chromosome"]==ii)
		for (i in 1:nrow(tmp)) {
			points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
		}
		for (i in 1:(nrow(tmp)-1)) {
			points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
		}
		abline(h=0, lwd=1)
		axis(2, at = c(-2,-1.5,-1,-.5,0,.5,1,1.5,2), labels=c("-2","-1.5","-1","-.5","0",".5","1","1.5","2"), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
		mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
		for (k in 1:8) {
			abline(h=(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), col="darkorange", lty=3)
			mtext(text=k, side=4, line=.5, at=(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), las=2, cex=.75, col="darkorange")
		}
		box(lwd=2)
		screen(zz[2])
		arg = copynumber:::getPlotParameters(type = "sample", nSeg = 10, cr = 3 * 3, sampleID = "dummy", plot.ideo = TRUE, xaxis = TRUE, assembly = "hg19")
		copynumber:::plotIdeogram(chrom=ii, TRUE, cyto.data = arg$assembly, cex = .75, unit = "bp")
		close.screen(all.screens=TRUE)
		dev.off()
	}


}
