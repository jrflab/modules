#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample", default = NA, type = 'character', help = "sample name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load("modules/copy_number/CytoBand.RData")
infile = paste0("qdnaseq/bed/", opt$sample, ".bed")
outfile = paste0("qdnaseq/copynumber/segmented/", opt$sample, ".RData")

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
data = read.table(file=infile, header=FALSE, sep="\t", skip=1, stringsAsFactors=FALSE)[,c(1,2,3,5),drop=FALSE]
colnames(data) = c("Chromosome", "Start", "End", "Log2Ratio")
segmented = pcf(data=winsorize(data=data[,c("Chromosome", "Start", "Log2Ratio"),drop=FALSE], method="mad", tau=2.5, k=25, verbose=FALSE), kmin = 250, gamma=150, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
colnames(segmented) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
save(data, segmented, file=outfile)
