#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("readr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_sets", default = NA, type = 'character', help = "sample sets file names"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = na.omit(unlist(strsplit(x=opt$sample_sets, split=" ", fixed=TRUE)))
list_of_dfs = list()
for (i in 1:length(sample_names)) {
	sample_vars = read_tsv(file=paste0("sufam/", sample_names[i], ".tsv"))
	col_names = colnames(sample_vars)
	sample_vars = as.data.frame(sample_vars)
	sample_vars[sample_vars=="" | sample_vars==" " | is.na(sample_vars)] = "NA"
	colnames(sample_vars) = col_names
	list_of_dfs[[i]] = sample_vars
}
names(list_of_dfs) = sample_names
write.xlsx(list_of_dfs, file="summary/sufam_summary.xlsx")
