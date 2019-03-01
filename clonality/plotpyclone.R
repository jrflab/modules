#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

optList = list(make_option("--sample_set", default = NULL, help = "sample set name"),
			   make_option("--normal_samples", default = NULL, help = "normal sample names"),
			   make_option("--burnin", default = NULL, help = "number of burnin mcmc"))
			   
parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_samples = unlist(strsplit(opt$sample_set, split="_", fixed=TRUE))
normal_sample = unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE))
normal_sample = tumor_samples[tumor_samples %in% normal_sample]
tumor_samples = tumor_samples[!(tumor_samples %in% normal_sample)]

ccf = list()
for (i in length(tumor_samples)) {
	data = read_tsv(file=paste0("pyclone/", opt$sample_set, "/trace/", tumor_samples[i], ".cellular_prevalence.tsv.bz2"))
	data = data[-(1:opt$burnin),,drop=FALSE]
	
	x = apply(data, 2, mean, na.rm=TRUE)

	if (dir.exists(paste0("pyclone/", opt$sample_set, "/plots"))) {
		dir.create(paste0("pyclone/", opt$sample_set, "/plots"))
	}
	pdf(paste0("pyclone/", opt$sample_set, "/plots/", tumor_samples[i], "_histogram_ccf.pdf"))
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	hist(x*100, col="grey90", border="grey80", axes=FALSE, main="", xlab="", ylab="", xlim=c(0,100))
	axis(1, at=NULL, cex.axis=1.5, padj=0.25)
    axis(2, at=NULL, cex.axis=1.5, las=1)
    mtext(side=1, text="Cancer cell fraction (%)", line=4, cex=1.5)
    mtext(side=2, text="Frequency", line=4, cex=1.5)
    dev.off()
    
    ccf[[i]] = x
}	

if (dir.exists(paste0("pyclone/", opt$sample_set, "/plots"))) {
	dir.create(paste0("pyclone/", opt$sample_set, "/plots"))
}
pdf(paste0("pyclone/", opt$sample_set, "/plots/2_by_2_scatter_plots.pdf"))
par(mar=c(6.1, 6.5, 4.1, 1.1))
for (i in 1:(length(ccf)-1)) {
	for (j in 2:(length(ccf))) {
		plot(ccf[[i]]*100, ccf[[j]]*100, pch=21, col="salmon3", bg="grey90", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", xlim=c(0,100), ylim=c(0,100))
    	axis(1, at=NULL, cex.axis=1.5, padj=0.25)
    	axis(2, at=NULL, cex.axis=1.5, las=1)
    	mtext(side=1, text=tumor_samples[i], line=4, cex=1.5)
    	mtext(side=2, text=tumor_samples[j], line=4, cex=1.5)
    	abline(h=10, lty=2, col="goldenrod3")
    	abline(h=20, lty=3, col="goldenrod3")
    	abline(h=90, lty=2, col="goldenrod3")
    	abline(h=80, lty=3, col="goldenrod3")
    	abline(v=10, lty=2, col="goldenrod3")
    	abline(v=20, lty=3, col="goldenrod3")
    	abline(v=90, lty=2, col="goldenrod3")
    	abline(v=80, lty=3, col="goldenrod3")
    	box(lwd=2)
    }
}
dev.off()
