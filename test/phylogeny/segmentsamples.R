#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
				  make_option("--normal_samples", default = NA, type = 'character', help = "normal samples"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_samples = na.omit(unlist(strsplit(opt$sample_set, split="_", fixed=TRUE)))
normal_samples = na.omit(unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE)))
normal_samples = normal_samples[normal_samples %in% all_samples]
tumor_samples = all_samples[!(all_samples %in% normal_samples)]

load(paste0("medicc/mad/", opt$sample_set, ".RData"))
# 
# 		multipcf(data, pos.unit = "bp", arms = NULL, Y = NULL, gamma = 40, 
#                normalize=TRUE, w=1, fast = TRUE, assembly = "hg19", digits = 4,
#                return.est = FALSE, save.res = FALSE, file.names = NULL, verbose 
#                = TRUE)

save(list=ls(all=TRUE), file=paste0("medicc/mpcf/", opt$sample_set, ".RData"))
