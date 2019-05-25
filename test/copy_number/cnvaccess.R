#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'plot_log2_' <- function(x, title = "")
{
   	par(mar=c(5, 5, 4, 2)+.1)
   	data("CytoBand")
   	end = NULL
   	for (i in 1:23) {
   		end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   	}
   	end = cumsum(end)
   	start = c(1, end[1:22]+1)
   	CytoBand = cbind(start, end)
   	index = NULL
   	for (i in 1:23) {
   		index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chromosome==i)))
   	}
	plot(index, x$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3", lty=3, lwd=.5)
	abline(h=0, col="red", lty=1, lwd=1)
	for (j in 1:23) {
		abline(v=CytoBand[j,"end"], col="goldenrod3", lty=3, lwd=.5)
	}
	axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = title, line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

if (as.numeric(opt$type)==1) {
	
	dataA = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".A.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataA[dataA[,"chromosome"]=="X", "chromosome"] = 23
	dataA[dataA[,"chromosome"]=="Y", "chromosome"] = 24
	dataA[,"chromosome"] = as.numeric(dataA[,"chromosome"])
	dataA = subset(dataA, dataA[,"chromosome"]<=23)
	
	tmp = dataA[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataA[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataB = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".B.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataB[dataB[,"chromosome"]=="X", "chromosome"] = 23
	dataB[dataB[,"chromosome"]=="Y", "chromosome"] = 24
	dataB[,"chromosome"] = as.numeric(dataB[,"chromosome"])
	dataB = subset(dataB, dataB[,"chromosome"]<=23)
	
	tmp = dataB[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataB[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataC = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".C.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataC[dataC[,"chromosome"]=="X", "chromosome"] = 23
	dataC[dataC[,"chromosome"]=="Y", "chromosome"] = 24
	dataC[,"chromosome"] = as.numeric(dataC[,"chromosome"])
	dataC = subset(dataC, dataC[,"chromosome"]<=23)
	
	tmp = dataC[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=1.5, k=25, verbose=FALSE, return.outliers=TRUE)
	dataC[tmp$wins.outliers[,3]!=0,"log2"] = NA

	data = rbind(dataA, dataB, dataC)
	index = order(data[,2])
	data = data[index,,drop=FALSE]
	index = order(data[,1])
	data = data[index,,drop=FALSE]
	
	pdf(file = paste0("cnvaccess/plot/log2/", opt$sample_name, ".pdf"), width=10, height=4.2)
	plot_log2_(x=data, title=opt$sample_name)
	dev.off()
	
}
