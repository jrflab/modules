#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--input", default = NA, type = 'character', help = "input file path"),
	       make_option("--output", default = NA, type = 'character', help = "output file path"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	smry = readr::read_tsv(file = opt$input, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::rename(`#CHROM` = CHROM,
			     POS = POS) %>%
	       dplyr::mutate(ID = ".",
			     QUAL = 100,
			     FILTER = "PASS",
			     INFO = ".") %>%
	       dplyr::select(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO)
	cat("##fileformat=VCFv4.2\n", file = opt$output, append = FALSE) 
	readr::write_tsv(smry, path = opt$output, na = "NA", append = TRUE, col_names = TRUE)
}
