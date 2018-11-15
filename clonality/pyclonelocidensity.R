#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--sample_name", default = NULL, help = "sample name"),
			   make_option("--burnin", default = NULL, help = "number of burnin iterations"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'transparent_rgb' <- function (col = "black", alpha = 85)
{
	tmp = c(col2rgb(col), alpha, 255)
	names(tmp) = c("red", "green", "blue", "alpha", "maxColorValue")
    out = do.call("rgb", as.list(tmp))
    return(invisible(out))
}

'post_density' <- function (x)
{
	y = density(x=x)
	return(invisible(y))
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
pdf(file=paste0("pyclone/", opt$sample_name, "/pyclone_loci_density.pdf"), height=5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
for (i in 1:length(feature_names)) {
	tmp = list()
	for (j in 1:length(ccf)) {
		tmp[[j]] = post_density(ccf[[j]][,i])
	}
	plot(0, 0, type="n", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0, max(unlist(lapply(tmp, function(x) { x$y })))))
	for (j in 1:length(tmp)) {
		points(tmp[[j]]$x, tmp[[j]]$y, col="black")
	}
    axis(1, at=NULL, cex.axis=1.5, padj=0.25)
    axis(2, at=NULL, cex.axis=1.5, las=1)
    mtext(side=1, text="CCF", line=4, cex=1.5)
    mtext(side=2, text="Density", line=4, cex=1.5)
    box(lwd=2)
}
dev.off()


