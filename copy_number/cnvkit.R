#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

'plot_log2_' <- function(x, y, title = "", alpha=NA, psi=NA)
{
   	par(mar=c(5, 5, 4, 2)+.1)
   	data("CytoBand")
   	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		y[y[,"Chromosome"]==j,"Start"] = y[y[,"Chromosome"]==j,"Start"] + start[j]
		y[y[,"Chromosome"]==j,"End"] = y[y[,"Chromosome"]==j,"End"] + start[j]
		x[x[,"chrom"]==j,"pos"] = x[x[,"chrom"]==j,"pos"] + start[j]
	}
	plot(x[,"pos"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey75", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
	for (j in 1:nrow(y)) {
 		lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=1.75, col="red")
 	}
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3", lty=3, lwd=.5)
	abline(h=0, col="red", lty=1, lwd=1)
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=.5)
	}
	abline(v=max(x[,"pos"]), col="goldenrod3", lty=3, lwd=.5)
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)	
	rect(xleft=1-1e10, xright=x[nrow(x),"pos"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = paste0(title, " | alpha = ", signif(alpha, 3), " | psi = ", signif(psi, 3)), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--type", default = NA, type = 'character', help = "type of analysis"),
				  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

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
	CN[,"Chromosome"] = gsub(pattern="chr", replacement="", x=CN[,"Chromosome"], fixed=TRUE)
	CN[CN[,"Chromosome"]=="X","Chromosome"] = 23
	CN[CN[,"Chromosome"]=="Y","Chromosome"] = 24
	CN[,"Chromosome"] = as.numeric(CN[,"Chromosome"])
	CN[CN[,"Log2Ratio"]<(-4) | CN[,"Log2Ratio"]>(4),"Log2Ratio"] = 0
	CN = subset(CN, CN[,"Chromosome"]<=23)
	tmp = pcf(data=winsorize(data=CN, method="mad", tau=2.5, k=10, verbose=FALSE), kmin = 10, gamma=40, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	save(CN, tmp, file=paste0("cnvkit/totalcopy/", opt$sample_name, ".RData"))
	tmp = prunesegments.cn(x=tmp, n=10)
	CN = winsorize(data=CN[,c("Chromosome","Position","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	pdf(file=paste0("cnvkit/segmented/", opt$sample_name, ".pdf"), width=10, height=4.25)
	file_names = dir(path="facets/cncf", pattern=opt$sample_name, full.names=TRUE)
	file_names = file_names[grep(".Rdata", file_names, fixed=TRUE)]
	if (length(file_names)==1) {
		load(file_names)
		alpha = fit$purity
		psi = fit$ploidy
	} else {
		alpha = NA
		psi = NA
	}
	plot_log2_(x=CN, y=tmp, title = opt$sample_name, alpha=alpha, psi=psi)
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
	file_names = dir(path="facets/cncf", pattern=opt$sample_name, full.names=TRUE)
	file_names = file_names[grep(".Rdata", file_names, fixed=TRUE)]
	if (length(file_names)==1) {
		load(file_names)
		alpha = ifelse(is.na(fit$purity), 1, fit$purity)
		psi = ifelse(is.na(fit$ploidy), 2, fit$ploid)
	} else {
		alpha = 1
		psi = 2
	}
	tmp = prunesegments.cn(x=tmp, n=10)
	qt = round((((2^(tmp[,"Log2Ratio"])) * (alpha*psi + 2*(1-alpha))) - 2*(1-alpha))/alpha)
	qt[is.na(qt)] = 2
	qt[is.infinite(qt)] = 2
	cat5 = rep(0, length(qt))
	if (round(psi)==1 | round(psi)==2) {
		cat5t = c(0, 1, 3, 7)
	} else if (round(psi)==3) {
		cat5t = c(0, 1, 4, 9)
	} else if (round(psi)==4) {
		cat5t = c(0, 1, 5, 10)
	} else if (round(psi)==5) {
		cat5t = c(0, 2, 6, 12)
	} else if (round(psi)>=6) {
		cat5t = c(0, 2, 7, 15)
	} else {
		cat5t = c(0, 1, 3, 7)
	}
	cat5[qt <= cat5t[2]] = -1
	cat5[qt <= cat5t[1]] = -2
	cat5[qt >= cat5t[3]] = 1
	cat5[qt >= cat5t[4]] = 2
	tmp = cbind(tmp, "Cat5"=cat5)
	save(CN, tmp, file=paste0("cnvkit/called/", opt$sample_name, ".RData"))
	
}

warnings()
