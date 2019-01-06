#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_names", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

sample_names = unlist(strsplit(x=opt$sample_names, split=" ", fixed=TRUE))
tsv = list()
for (i in 1:length(sample_names)) {
	tsv[[i]] = read.csv(file=paste0("cravat/", sample_names[i], ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
tsv = do.call(rbind, tsv)
write.table(tsv, file="summary/tsv/cravat_summary.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
