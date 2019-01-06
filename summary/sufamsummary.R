#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("readr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

file_names = unlist(strsplit(x=opt$in_file, split=" ", fixed=TRUE))
sample_names = gsub(".tsv", "", x=gsub("sufam/", "", file_names, fixed=TRUE), fixed=TRUE)
list_of_dfs = list()
for (i in 1:length(file_names)) {
	sample_vars = read_tsv(file=file_names[i])
	col_names = colnames(sample_vars)
	sample_vars = as.data.frame(sample_vars)
	colnames(sample_vars) = col_names
	list_of_dfs[[i]] = sample_vars
}
names(list_of_dfs) = sample_names
write.xlsx(list_of_dfs, file="summary/sufam_summary.xlsx")
