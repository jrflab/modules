#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--type", default = NA, type = 'character', help = "type of analysis"),
				  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options
load("modules/copy_number/CytoBand.RData")

if (opt$type=="total-copy") {
	
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
		return(x)
	}

	data = read.csv(file=paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	CN = data[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
	CN[CN[,"Chromosome"]=="X","Chromosome"] = 23
	CN[CN[,"Chromosome"]=="Y","Chromosome"] = 24
	CN[,"Chromosome"] = as.numeric(CN[,"Chromosome"])
	CN[CN[,"Log2Ratio"]<(-4) | CN[,"Log2Ratio"]>(4),"Log2Ratio"] = 0
	CN = subset(CN, CN[,"Chromosome"]<=23)
	tmp = pcf(data=winsorize(data=CN, method="mad", tau=2.5, k=25, verbose=FALSE), kmin = 50, gamma=70, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	save(CN, tmp, file=paste0("cnvkit/totalcopy/", opt$sample_name, ".RData"))
	tmp = prunesegments.cn(x=tmp, n=10)
	CN = winsorize(data=CN[,c("Chromosome","Position","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		tmp[tmp[,"Chromosome"]==j,"Start"] = tmp[tmp[,"Chromosome"]==j,"Start"] + start[j]
		tmp[tmp[,"Chromosome"]==j,"End"] = tmp[tmp[,"Chromosome"]==j,"End"] + start[j]
		CN[CN[,"chrom"]==j,"pos"] = CN[CN[,"chrom"]==j,"pos"] + start[j]
	}
	col = "grey70"
	pdf(file=paste0("cnvkit/segmented/", opt$sample_name, ".pdf"), height=7*10/7*2/3, width=7*20/7)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(CN[,"pos"], CN[,"Log2Ratio"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
 	for (j in 1:nrow(tmp)) {
 		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
 	}
 	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(CN[,"pos"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	dev.off()
	

} else if (opt$type=="call-cna") {
	
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
		return(x)
	}
	load(paste0("cnvkit/totalcopy/", opt$sample_name, ".RData"))
	tmp = prunesegments.cn(x=tmp, n=10)
	cat5t = c(-0.4, -0.19, 0.15, 0.58)
	cat5 = rep(0, nrow(tmp))
	cat5[tmp[,"Log2Ratio"] < cat5t[2]] = -1
	cat5[tmp[,"Log2Ratio"] < cat5t[1]] = -2
	cat5[tmp[,"Log2Ratio"] > cat5t[3]] = 1
	cat5[tmp[,"Log2Ratio"] > cat5t[4]] = 2
	tmp = cbind(tmp, "Cat5"=cat5)
	save(CN, tmp, file=paste0("cnvkit/calls/", opt$sample_name, ".RData"))
	
}

warnings()
