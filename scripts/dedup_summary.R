#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".dedup_metrics.txt"),
					 skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="metrics/dedup_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

}
