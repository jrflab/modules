#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "type of analysis"),
		  make_option("--sample_names", default = NA, type = 'character', help = "sample name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option)==1) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	smry = list()
	for (i in 1:length(sample_names)) {
		smry[[i]] = readr::read_tsv(file = paste0("mimsi/", sample_names[i], "/", sample_names[i], ".txt"),
					    col_names = TRUE, col_types = cols(.default = col_character())) %>%
			    readr::type_convert()
	}
	smry = do.call(rbind, smry)
	write_tsv(smry, path="mimsi/summary.txt", append = FALSE, col_names = TRUE)
	
}
