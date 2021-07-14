#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file_name", default = NA, type = 'character', help = "sample names set"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

genotype = readr::read_tsv(file = opt$file_name, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(Chrom_N = gsub(pattern = "chr", replacement = "" x = Chrom, fixed = TRUE)) %>%
	   readr::type_convert() %>%
	   dplyr::arrange(Chrom_N, Pos) %>%
	   dplyr::select(-Chrom_N)

write_tsv(genotype, file = gsub(pattern = ".txt", replacement = ".tsv", x = opt$file_name), append = FALSE, col_names = TRUE)
