#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("Starfish"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_name = as.character(opt$sample_name)
	bed = readr::read_tsv(file = as.character(opt$input_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::rename(chrom1 = X1,
			    start1 = X2,
			    end1 = X3,
			    chrom2 = X4,
			    start2 = X5,
			    end2 = X6,
			    sv_id = X7,
			    pe_support = X8,
			    strand1 = X9,
			    strand2 = X10,
			    svclass = X11) %>%
	      dplyr::select(chrom1, pos1 = start1, chrom2, pos2 = start2, strand1, strand2, svtype = svclass) %>%
	      dplyr::mutate(svtype = case_when(
		      svtype == "INV" & strand1 == "+" & strand2 == "+" ~ "h2hINV",
		      svtype == "INV" & strand1 == "-" & strand2 == "-" ~ "t2tINV",
		      TRUE ~ svtype
	      )) %>%
	      dplyr::mutate(sample = sample_name)
	readr::write_tsv(x = bed, file = as.character(opt$output_file), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==2) {
	sample_name = as.character(opt$sample_name)
	data = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       dplyr::select(chromosome = chrom,
			     start = loc.start,
			     end = loc.end,
			     total_cn = tcn.em) %>%
	       dplyr::mutate(sample = sample_name)
	readr::write_tsv(x = data, file = as.character(opt$output_file), col_names = TRUE, append = FALSE)
	
}

