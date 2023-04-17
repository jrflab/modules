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
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".idx_stats.txt"),
					 col_names = FALSE, col_types = cols(.default = col_character()))[-85,,drop=FALSE] %>%
			 readr::type_convert() %>%
			 dplyr::select(CHROMOSOME = X1,
				       LENGTH = X2,
				       ALIGNED_READS = X3) %>%
			 dplyr::mutate(CHROMOSOME = gsub(pattern=" length=", replacement="", x=CHROMOSOME),
				       ALIGNED_READS = gsub(pattern="Aligned= ", replacement="", x=ALIGNED_READS),
				       SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/idx_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==2) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".aln_metrics.txt"),
					 skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::select(-SAMPLE, -READ_GROUP) %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/aln_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==3) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".insert_metrics.txt"),
					 skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::select(-SAMPLE, -READ_GROUP) %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/insert_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==4) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".oxog_metrics.txt"),
					 skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::rename(SAMPLE_NAME = SAMPLE_ALIAS)
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/oxog_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==5) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".hs_metrics.txt"),
					 skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/hs_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==6) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".gc_metrics.txt"),
					 skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/gc_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==7) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		x[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".gc_bias.txt"),
					 skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x = do.call(rbind, x)
	write_tsv(x, path="summary/gc_summary.txt", na = "NA", append = FALSE, col_names = TRUE)

}
