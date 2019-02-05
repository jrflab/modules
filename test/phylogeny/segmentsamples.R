#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

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
save(list=ls(all=TRUE), file=paste0("medicc/mad/", opt$sample_set, ".RData"))
