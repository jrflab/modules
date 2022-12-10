#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "Which option?"),
		  make_option("--tumor_name", default = NA, type = 'character', help = "Tumor sample name"),
		  make_option("--normal_name", default = NA, type = 'character', help = "Normal sample name"),
		  make_option("--file_in", default = NA, type = 'character', help = "Input file name including path"),
		  make_option("--file_out", default = NA, type = 'character', help = "Output file name including path"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option) == 1) {
	load(as.character(opt$file_in))
	cn_df = out2$jointseg %>%
		dplyr::as_tibble() %>%
		dplyr::filter(het == 1) %>%
		dplyr::select(Chromosome = chrom,
			      Position = maploc,
			      Log2_Ratio = cnlr,
			      B_Allele_Frequency = vafT)
	readr::write_tsv(x = cn_df, file = as.character(opt$_file_out), col_names = TRUE, append = FALSE)

}