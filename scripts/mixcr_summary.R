#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

sample_names = unlist(strsplit(x = as.character(opt$sample_names), split = " ", fixed = TRUE))
clones = list()
for (i in 1:length(sample_names)) {
	file_name = paste0("mixcr/", sample_names[i], "/clones.tsv")
	if (!file.exists(file_name)) {
		file_names = dir(path = paste0("mixcr/", sample_names[i]), pattern = "clones_")
		df = list()
		for (j in 1:length(file_names)) {
			df[[j]] = readr::read_tsv(file = paste0("mixcr/", sample_names[i], "/", file_names[j]), col_names = TRUE, col_types = cols(.default = col_character())) %>%
					  readr::type_convert() %>%
					  dplyr::mutate(sample_name = sample_names[i]) %>%
					  dplyr::mutate(clone = gsub(pattern = "clones_|.tsv", replacement = "", x = file_names[j], fixed = FALSE, perl = TRUE))
		}
		clones[[i]] = do.call(rbind, df)
	}
}
clones = do.call(rbind, clones)
readr::write_tsv(x = clones, path = "mixcr/summary.txt", append = FALSE, col_names = TRUE)