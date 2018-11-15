#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--sample_name", default = NULL, help = "sample name"),
			   make_option("--burnin", default = NULL, help = "number of burnin iterations"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'post_density' <- function (x)
{
	y = density(x=x)
	return(invisible(y))
}

'hex_cols' <- function(x)
{
	cols = c("#4865B1", "#FFA500", "#B22034", "#E9E0Ba", "#D5D5D5", "#000000", "#DC0073", "#00A1E4")
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
pdf(file=paste0("pyclone/", opt$sample_name, "/plots/by_loci_density.pdf"), width=8)
par(mar=c(6.1, 6.5, 4.1, 1.1))
for (i in 1:length(feature_names)) {
	tmp = list()
	for (j in 1:length(ccf)) {
		tmp[[j]] = post_density(ccf[[j]][,i])
	}
	plot(0, 0, type="n", axes=FALSE, frame.plot=FALSE, main=sub("_", " ", feature_names[i], fixed=TRUE), xlab="", ylab="", xlim=c(0,1), ylim=c(0,1.1), cex.main=2)
	for (j in 1:length(tmp)) {
		index = tmp[[j]]$x>1 | tmp[[j]]$x<0
		points(tmp[[j]]$x[!index], ((tmp[[j]]$y-min(tmp[[j]]$y))/(max(tmp[[j]]$y)-min(tmp[[j]]$y)))[!index], type="l", lwd=3, col=hex_cols(j))
		points(rep(mean(ccf[[j]][,i]),2), c(-1,1), type="l", col=hex_cols(j), lty=2, lwd=1)
	}
    axis(1, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, padj=0.25, lwd = 1.25, lwd.ticks = 1.15)
    axis(2, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, las=1, lwd = 1.25, lwd.ticks = 1.15)
    mtext(side=1, text="CCF", line=4, cex=1.5)
    mtext(side=2, text="Density", line=4, cex=1.5)
}
dev.off()
