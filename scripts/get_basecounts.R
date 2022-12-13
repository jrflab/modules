#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "Which option?"),
		  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option) == 1) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_name), split = " ", fixed = TRUE))
	data = list()
	for (i in 1:length(sample_names)) {
		data[[i]] = readr::read_tsv(file = paste0("gbc/", sample_names[i], ".txt.gz"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			    readr::type_convert() %>%
			    dplyr::mutate(sample_name = sample_names[i])
	}
	data = do.call(bind_rows, data)
	readr::write_tsv(x = data, path = "gbc/summary.txt", append = FALSE, col_names = TRUE)
}

