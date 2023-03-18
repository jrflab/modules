#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--input_file", default = NA, type = 'character', help = "Input VCF file"),
               make_option("--output_file", default = NA, type = 'character', help = "Output VCF file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

vcf = readr::read_tsv(file = as.character(opt$input_file), comment = "#", col_names = FALSE, col_types = cols(.default = col_character())) %>%
      readr::type_convert() %>%
      dplyr::filter(!grepl("SUPP_VEC=110", X8, fixed = TRUE)) %>%
      dplyr::mutate(X3 = X12) %>%
      dplyr::mutate(X3 = unlist(lapply(X3, function(x) { unlist(strsplit(x, split = ":", fixed = TRUE))[8] }))) %>%
      dplyr::mutate(X3 = gsub(pattern = "_", replacement = ":", x = X3, fixed = TRUE)) %>%
      dplyr::mutate(X5 = case_when(
	      grepl("DUP", X3, fixed = TRUE) ~ "<DUP:TANDEM>",
	      grepl("DEL", X3, fixed = TRUE) ~ "<DEL>",
	      grepl("INV", X3, fixed = TRUE) ~ "<INV>",
	      TRUE ~ X5
      )) %>%
      dplyr::mutate(X8 = case_when(
	      grepl("DUP", X3, fixed = TRUE) ~ gsub("SVTYPE=INV", "SVTYPE=DUP", X8),
	      grepl("DEL", X3, fixed = TRUE) ~ gsub("SVTYPE=INV", "SVTYPE=DEL", X8),
	      TRUE ~ X8
      )) %>%
      dplyr::rename(`#CHROM` = X1,
		    POS = X2,
		    ID = X3,
		    REF = X4,
		    ALT = X5,
		    QUAL = X6,
		    FILTER = X7,
		    INFO = X8,
		    FORMAT = X9,
		    SVABA = X10,
		    GRIDSS = X11,
		    MANTA = X12)

readr::write_tsv(x = vcf, path = as.character(opt$output_file), append = TRUE, col_names = TRUE)


