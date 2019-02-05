#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
				  make_option("--normal_samples", default = NA, type = 'character', help = "sample names set"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_samples = na.omit(unlist(strsplit(opt$sample_set, split="_", fixed=TRUE)))
normal_samples = na.omit(unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE)))
normal_samples = normal_samples[normal_samples %in% all_samples]
tumor_samples = all_samples[!(all_samples %in% normal_samples)]

CN = list()
for (i in 1:length(tumor_samples)) {
	load(paste0("facets/cncf/", tumor_samples[i], "_", normal_samples, ".Rdata"))
	CN[[i]] = out2$jointseg[,c("chrom", "maploc", "cnlr"),drop=FALSE]
	colnames(CN[[i]]) = c("Chromosome", "Position", "Log2Ratio")
}
save(CN, file=paste0("medicc/ascat/", opt$sample_set, ".RData"))
