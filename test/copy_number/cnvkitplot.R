#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

'MAPD' <- function(x) {
	i = 1
	ii = length(x)-1
	j = 2
	jj = length(x)
	y = median(abs(x[i:ii] - x[j:jj]))
	return(invisible(y))
}

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
outfile_on_target_A = paste0("cnvkit/plot/", tumor_sample, ".A.ontarget.pdf")
outfile_on_target_B = paste0("cnvkit/plot/", tumor_sample, ".B.ontarget.pdf")
outfile_on_target_AB = paste0("cnvkit/plot/", tumor_sample, ".AB.ontarget.pdf")
outfile_off_target = paste0("cnvkit/plot/", tumor_sample, ".offtarget.pdf")

save = list()

data = list()
for (i in 1:length(normal_samples)) {
	file = paste0("cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".A.cnr")
	data[[i]] = read.table(file=file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
}
m0 = vector(mode="numeric", length(normal_samples))
for (i in 1:length(normal_samples)) {
	index = data[[i]][,"gene"] == "-"
	m0[i] = MAPD(data[[i]][index,"log2"])
}
index = data[[which.min(m0)]][,"gene"]=="-"
data = data[[which.min(m0)]][index,,drop=FALSE]
save[[1]] = data
if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target_A))
} else {
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<23)
	col = rep("#9F6986", nrow(data))
	col[(data[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_on_target_A, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"log2"], type="p", pch=19, cex=.25, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
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
	box(lwd=2.5)
	dev.off()
}

data = list()
for (i in 1:length(normal_samples)) {
	file = paste0("cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".B.cnr")
	data[[i]] = read.table(file=file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
}
m0 = vector(mode="numeric", length(normal_samples))
for (i in 1:length(normal_samples)) {
	index = data[[i]][,"gene"] == "-"
	m0[i] = MAPD(data[[i]][index,"log2"])
}
index = data[[which.min(m0)]][,"gene"]=="-"
data = data[[which.min(m0)]][index,,drop=FALSE]
save[[2]] = data
if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target_B))
} else {
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<23)
	col = rep("#9F6986", nrow(data))
	col[(data[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_on_target_B, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"log2"], type="p", pch=19, cex=.25, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
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
	box(lwd=2.5)
	dev.off()
}

data = list()
for (i in 1:length(normal_samples)) {
	file = paste0("cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".A.cnr")
	data[[i]] = read.table(file=file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
}
m0 = vector(mode="numeric", length(normal_samples))
for (i in 1:length(normal_samples)) {
	index = data[[i]][,"gene"] != "-"
	m0[i] = MAPD(data[[i]][index,"log2"])
}
index = data[[which.min(m0)]][,"gene"]!="-"
data = data[[which.min(m0)]][index,,drop=FALSE]
save[[3]] = data
if (nrow(data)==0) {
	system(paste0("touch ", outfile_off_target))
} else {	
	data = subset(data, data$gene!="-")
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<23)
	tmp = data[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(tmp) = c("Chromosome", "Position", "Log2Ratio")
	tmp = winsorize(data=tmp, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
	data[tmp$wins.outliers[,3]!=0,"log2"] = NA
	col = rep("#9F6986", nrow(data))
	col[(data[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_off_target, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"log2"], type="p", pch=".", cex=2, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
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
	box(lwd=2.5)
	dev.off()
}

data = rbind(save[[1]], save[[2]])
if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target_AB))
} else {
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<23)
	index = order(data[,"start"])
	data = data[index,,drop=FALSE]
	index = order(data[,"chromosome"])
	data = data[index,,drop=FALSE]
	data[data[,"log2"]<(-2),"log2"] = 0
	tmp = save[[3]][,c("chromosome", "start", "end", "log2"),drop=FALSE]
	tmp[tmp[,"chromosome"]=="X", "chromosome"] = 23
	tmp[tmp[,"chromosome"]=="Y", "chromosome"] = 24
	tmp[,"chromosome"] = as.numeric(tmp[,"chromosome"])
	tmp = subset(tmp, tmp[,"chromosome"]<23)
	index = order(tmp[,"start"])
	tmp = tmp[index,,drop=FALSE]
	index = order(tmp[,"chromosome"])
	tmp = tmp[index,,drop=FALSE]
	tmp[,"log2"] = NA
	for (j in 1:22) {
		ind0 = which(tmp[,"chromosome"]==j)
		ind1 = which(data[,"chromosome"]==j)
		set.seed(0)
		indx = sort(sample(x=ind0, size=length(ind1)))
		tmp[indx,"log2"] = data[ind1,"log2"]
	}
	data = tmp
	col = rep("#9F6986", nrow(data))
	col[(data[,"chromosome"]%%2)==1] = "#CECAC5"
	pdf(file=outfile_on_target_AB, width=14, height=5)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(data[,"log2"], type="p", pch=19, cex=.25, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
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
	box(lwd=2.5)
	dev.off()
}
