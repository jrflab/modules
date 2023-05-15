#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "type of analysis"),
		  make_option("--sample_pairs", default = NA, type = 'character', help = "sample pairs"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option) == 1) {
	sample_names = unlist(strsplit(as.character(opt$sample_pairs), split = " ", fixed = TRUE))
	CN = list()
	for (i in 1:length(sample_names)) {
		CN[[i]] = readr::read_tsv(file = paste0("facets_suite/", sample_names[i], "/", sample_names[i], ".gene_level.txt"),
					  col_names = TRUE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert()
	}
	CN = do.call(rbind, CN)
	readr::write_tsv(x = CN, path = "facets_suite/summary.txt", col_names = TRUE, append = FALSE)

}
