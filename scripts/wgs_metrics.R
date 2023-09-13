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
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".idx_stats.txt"),
					       col_names = FALSE, col_types = cols(.default = col_character()))[-85,,drop=FALSE] %>%
			       readr::type_convert() %>%
			       dplyr::select(CHROMOSOME = X1,
					     LENGTH = X2,
					     ALIGNED_READS = X3) %>%
			       dplyr::mutate(CHROMOSOME = gsub(pattern=" length=", replacement="", x=CHROMOSOME),
					     ALIGNED_READS = gsub(pattern="Aligned= ", replacement="", x=ALIGNED_READS),
					     SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/idx_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==2) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".aln_metrics.txt"),
					       skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::select(-SAMPLE, -READ_GROUP) %>%
			       dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/aln_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==3) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".insert_metrics.txt"),
					       skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::select(-SAMPLE, -READ_GROUP) %>%
			       dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/insert_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==4) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".oxog_metrics.txt"),
					  skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::rename(SAMPLE_NAME = SAMPLE_ALIAS)
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/oxog_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==5) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".gc_metrics.txt"),
					       skip = 6, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/gc_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==6) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".wgs_metrics.txt"),
					       skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/wgs_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==7) {
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	metrics = list()
	for (i in 1:length(sample_names)) {
		metrics[[i]] = readr::read_tsv(file = paste0("metrics/", sample_names[i], ".duplicate_metrics.txt"),
					       skip = 6, n_max = 1, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	metrics = do.call(rbind, metrics)
	write_tsv(metrics, path="summary/duplicate_metrics.txt", na = "NA", append = FALSE, col_names = TRUE)

}
