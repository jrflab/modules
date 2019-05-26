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
	tmp = tmp[index,,drop=FASE]
	index = order(tmp[,"Chromosome"])
	tmp = tmp[index,,drop=FASE]

	sample_name = gsub("_A.gz", "", x=strsplit(opt$pool_A, split="/", fixed=TRUE)[[1]][3], fixed=TRUE)
	write.table(tmp, file=paste0("fastaccess/cnr/", sample_name, ".txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
}

