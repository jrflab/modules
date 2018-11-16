#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("MASS"))

optList = list(make_option("--sample_name", default = NULL, help = "sample name"),
			   make_option("--burnin", default = NULL, help = "number of burnin iterations"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'post_density' <- function (x)
{
	y = density(x=x, adjust=2)
	return(invisible(y))
}

'hex_cols' <- function(x)
{
	x = x%%8
	if (x==0) {
		x = 8
	}
	cols = c("#4865B1", "#FFA500", "#B22034", "#E9E0BA", "#D5D5D5", "#000000", "#DC0073", "#00A1E4")
	return(cols[x])
}

file_names = dir(path=paste0("pyclone/", opt$sample_name, "/trace"), pattern="cellular_prevalence.tsv.bz2", full.names=TRUE)
ccf = list()
for (i in 1:length(file_names)) {
	ccf[[i]] = read.csv(file=file_names[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
feature_names = colnames(ccf[[1]])
for (i in 1:length(file_names)) {
	ccf[[i]] = ccf[[i]][,feature_names,drop=FALSE]
}
for (i in 1:length(file_names)) {
	ccf[[i]] = ccf[[i]][-(1:opt$burnin),,drop=FALSE]
}
pdf(file=paste0("pyclone/", opt$sample_name, "/plots/all_loci_scatter.pdf"))
par(mar=c(6.1, 6.5, 4.1, 1.1))
zz = matrix(NA, nrow=length(feature_names), ncol=length(ccf))
for (i in 1:length(ccf)) {
	z = vector(mode="numeric", length=length(feature_names))
	for (j in 1:length(feature_names)) {
		x = ccf[[i]][,j]
		y = Mclust(x, G=2)
		a = x[y$classification==1]
		b = x[y$classification==2]
		if (length(a)>length(b)) {
			pd = post_density(a)
			idx = which.max(pd$y)
			z[j] = pd$x[idx]
		} else {
			pd = post_density(b)
			idx = which.max(pd$y)
			z[j] = pd$x[idx]
		}
	}
	zz[,i] = z
}
for (i in 1:(ncol(zz)-1)) {
	for (j in (i+1):ncol(zz)) {
		plot(0, 0, type="n", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
		ind = is.na(zz[,i]) | is.na(zz[,j])
		x = jitter(zz[!ind,i], amount=.01)
		x[x>1] = 1
		x[x<0] = 0
		y = jitter(zz[!ind,j], amount=.01)
		y[y>1] = 1
		y[y<0] = 0
		contour(kde2d(x, y, n=100, lims = c(c(-.01,1.01),c(-.01,1.01))), drawlabels=FALSE, nlevels=15, add=TRUE, lwd=.5, col=hex_cols(1))
		points(x, y, type="p", pch=1, col=hex_cols(1))
	    axis(1, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, padj=0.25, lwd = 1.25, lwd.ticks = 1.15)
	    axis(2, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, las=1, lwd = 1.25, lwd.ticks = 1.15)
	    points(c(.1,.1), c(-.1,1), type="l", col="orange", lty=3)
	    points(c(.9,.9), c(-.1,1), type="l", col="orange", lty=3)
	    points(c(-.1,1), c(.1,.1), type="l", col="orange", lty=3)
	    points(c(-.1,1), c(.9,.9), type="l", col="orange", lty=3)
	    mtext(side=1, text=gsub(pattern="trace/", replacement="", x=gsub(pattern=paste0("pyclone/", opt$sample_name, "/"), replacement="", x=gsub(pattern=".cellular_prevalence.tsv.bz2", replacement="", x=file_names[i], fixed=TRUE), fixed=TRUE), fixed=TRUE), line=4, cex=1.5)
	    mtext(side=2, text=gsub(pattern="trace/", replacement="", x=gsub(pattern=paste0("pyclone/", opt$sample_name, "/"), replacement="", x=gsub(pattern=".cellular_prevalence.tsv.bz2", replacement="", x=file_names[j], fixed=TRUE), fixed=TRUE), fixed=TRUE), line=4, cex=1.5)
	}
	    
}
dev.off()
