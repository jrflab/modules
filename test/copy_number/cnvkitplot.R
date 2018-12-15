#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--tumor", default = NA, type = 'character', help = "tumor sample"),
				 make_option("--normals", default = NA, type = 'character', help = "normal samples"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_sample = opt$tumor
normal_samples = unlist(strsplit(x=opt$normals, " ", fixed=TRUE))
outfile_on_target = paste0("cnvkit/plot/", tumor_sample, ".ontarget.pdf")
outfile_off_target = paste0("cnvkit/plot/", tumor_sample, ".offtarget.pdf")

data = list()
for (i in 1:length(normal_samples)) {
	file = paste0("cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".cnr")
	data[[i]] = read.table(file=file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
}
mad0 = mad1 = vector(mode="numeric", length(normal_samples))
for (i in 1:length(normal_samples)) {
	index = data[[i]][,"gene"] == "-" & data[[i]][,"depth"]<1.5
	mad0[i] = mad(data[[i]][index,"log2"])
	index = data[[i]][,"gene"] != "-"
	mad1[i] = mad(data[[i]][index,"log2"])
}
index = data[[which.min(mad0)]][,"gene"]=="-"
data0 = data[[which.min(mad0)]][index,,drop=FALSE]
index = data[[which.min(mad1)]][,"gene"]!="-"
data1 = data[[which.min(mad1)]][index,,drop=FALSE]
data = rbind(data0, data1)
if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target))
	system(paste0("touch ", outfile_off_target))
} else {
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<=23)
	
	ontarget = subset(data, data$gene=="-" & data$depth<1.5)
	tmp = ontarget[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(tmp) = c("Chromosome", "Position", "Log2Ratio")
	tmp = winsorize(data=tmp, tau=3.5, k=5, verbose=FALSE, return.outliers=TRUE)
	ontarget[tmp$wins.outliers[,3]!=0,"log2"] = NA
	col = rep("#9F6986", nrow(ontarget))
	col[(ontarget[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_on_target, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(ontarget[,"log2"], type="p", pch=19, cex=.5, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:max(ontarget[,"chromosome"])) {
		v = min(which(ontarget[,"chromosome"]==j))
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(nrow(ontarget)), col="goldenrod3")
	start = NULL
	end = NULL
	for (j in 1:max(ontarget[,"chromosome"])) {
		start[j] = min(which(ontarget[,"chromosome"]==j))
		end[j] = max(which(ontarget[,"chromosome"]==j))
	}
	labels = 1:max(ontarget[,"chromosome"])
	labels[labels==23] = "X"
	axis(1, at = .5*(start+end), labels=labels, cex.axis = 0.85, las = 1)
	box(lwd=2.5)
	dev.off()
	
	offtarget = subset(data, data$gene!="-" & depth<1.5)
	tmp = offtarget[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(tmp) = c("Chromosome", "Position", "Log2Ratio")
	tmp = winsorize(data=tmp, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
	offtarget[tmp$wins.outliers[,3]!=0,"log2"] = NA
	col = rep("#9F6986", nrow(offtarget))
	col[(offtarget[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_off_target, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(offtarget[,"log2"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:max(offtarget[,"chromosome"])) {
		v = min(which(offtarget[,"chromosome"]==j))
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(nrow(offtarget)), col="goldenrod3")
	start = NULL
	end = NULL
	for (j in 1:max(offtarget[,"chromosome"])) {
		start[j] = min(which(offtarget[,"chromosome"]==j))
		end[j] = max(which(offtarget[,"chromosome"]==j))
	}
	labels = 1:max(offtarget[,"chromosome"])
	labels[labels==23] = "X"
	axis(1, at = .5*(start+end), labels=labels, cex.axis = 0.85, las = 1)
	box(lwd=2.5)
	dev.off()
}
