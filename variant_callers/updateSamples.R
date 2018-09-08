#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--patient", default = NA, type = 'character', help = "type of analysis"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_vars = read.csv(file=paste0("sufam/", opt$patient, ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
write.csv(all_vars, file=paste0("sufam/", opt$patient, ".tsv"))
