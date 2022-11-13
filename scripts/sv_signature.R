#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_names", default = NA, type = 'character', help = "sample names"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_names), split = " ", fixed=TRUE))
	feature_counts = list()
	for (i in 1:length(sample_names)) {
		feature_counts[[i]] = readr::read_tsv(file = paste0("sv_signature/", sample_names[i], "/", sample_names[i], ".merged.txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
				      readr::type_convert() %>%
				      dplyr::rename(sv_class = X1,
					            sv_count = manual_sv_type) %>%
				      dplyr::mutate(sample_name = sample_names[i])
	}
	feature_counts = do.call(bind_rows, feature_counts)
	write_tsv(x = feature_counts, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
}
