#!/usr/bin/env Rscript

set.seed(0)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("facets"))
suppressPackageStartupMessages(library("pctGCdata"))
suppressPackageStartupMessages(library("foreach"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--option", default = NA, type = 'numeric', help = "Type of analysis"),
                make_option("--sample_name", default = NA, type = 'character', help = "Sample name"),
                make_option("--pool_A", default = NA, type = 'character', help = "Pool A snp pileup"),
                make_option("--pool_B", default = NA, type = 'character', help = "Pool B snp pileup")
                )
parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option)==1) {
	snpmat_A = readSnpMatrix(filename=opt$pool_A)
	snpmat_B = readSnpMatrix(filename=opt$pool_B)
	proc_A = preProcSample(rcmat = snpmat_A,
						   ndepth = 50,
						   het.thresh = 0.1,
						   snp.nbhd = 10,
						   cval = 25,
						   deltaCN = 0,
						   gbuild = "hg19",
						   hetscale = TRUE,
						   unmatched = TRUE,
						   ndepthmax = 100000)
	proc_B = preProcSample(rcmat = snpmat_B,
						   ndepth = 25,
						   het.thresh = 0.1,
						   snp.nbhd = 10,
						   cval = 25,
						   deltaCN = 0,
						   gbuild = "hg19",
						   hetscale = TRUE,
						   unmatched = TRUE,
						   ndepthmax = 100000)
	tmp = rbind(proc_A$jointseg, proc_B$jointseg)[,c("chrom", "maploc", "cnlr", "vafT"),drop=FALSE]
	colnames(tmp) = c("Chromosome", "Position", "Log2Ratio", "BAF")
	index = order(tmp[,"Position"])
	tmp = tmp[index,,drop=FALSE]
	index = order(tmp[,"Chromosome"])
	tmp = tmp[index,,drop=FALSE]

	sample_name = gsub("_A.gz", "", x=strsplit(opt$pool_A, split="/", fixed=TRUE)[[1]][3], fixed=TRUE)
	write.table(tmp, file=paste0("fastaccess/cnr/", sample_name, ".txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
	
} else if (as.numeric(opt$option)==2) {
	
	'plot_log2_' <- function(x, title = "")
	{
		load("modules/copy_number/CytoBand.RData")
		par(mar=c(5, 5, 4, 2)+.1)
		end = NULL
		for (i in 1:23) {
			end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
		}
		end = cumsum(end)
		start = c(1, end[1:22]+1)
		CytoBand = cbind(start, end)
		index = NULL
		for (i in 1:23) {
			index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$Chromosome==i)))
		}
		plot(index, x$Log2Ratio, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		abline(v=1, col="goldenrod3", lty=3, lwd=.5)
		for (j in 1:23) {
			abline(v=CytoBand[j,"end"], col="goldenrod3", lty=3, lwd=.5)
		}
		axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
		rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
		title(main = title, line=-1, cex.main=.75, font.main=1)
		box(lwd=1.5)
	}
	
	print(opt$sample_name)
	data = read.csv(file=paste0("fastaccess/cnr/", opt$sample_name, ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	pdf(file=paste0("fastaccess/plots/log2/", opt$sample_name, ".pdf"), width=10, height=4.25)
	plot_log2_(x=data, title = opt$sample_name)
	dev.off()

}

