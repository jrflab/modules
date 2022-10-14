#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
	       make_option("--bar_code", default = NA, type = 'character', help = "sample name"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_name = opt$sample_name
	bar_code = opt$bar_code
	pile_up = readr::read_tsv(file = paste0("gbc/", sample_name, "/", bar_code, ".txt.gz"),
				  col_names = TRUE,
				  col_types = cols(.default = col_character())) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(A = A + a,
			        T = T + t,
			        G = G + g,
			        C = C + c) %>%
		  dplyr::rename(chromosome = Chrom,
			        position = Pos,
			        reference_allele = Ref,
			        total_depth = TOTAL_depth) %>%
		  dplyr::select(chromosome, position, reference_allele, total_depth, A, T, G, C) %>%
		  reshape2::melt(id.vars = c("chromosome", "position", "reference_allele", "total_depth"),
				 measure.vars = c("A", "T", "G", "C"),
				 variable.name = "alternate_allele",
				 value.name = "alternate_depth") %>%
		  dplyr::filter(reference_allele != alternate_allele) %>%
		  dplyr::group_by(chromosome, position) %>%
		  dplyr::summarize(reference_allele = unique(reference_allele),
				   total_depth = unique(total_depth),
				   alternate_depth = sum(alternate_depth)) %>%
		  dplyr::ungroup()
	write_tsv(pile_up, path = paste0("summary/", sample_name, "/sum_alt/", bar_code, ".txt"), na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==2) {
	sample_name = opt$sample_name
	pile_up = readr::read_tsv(file = paste0("gbc/", sample_name, "/", bar_code, ".txt.gz"),
				  col_names = TRUE,
				  col_types = cols(.default = col_character())) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(alternate_depth = INS + DEL) %>%
		  dplyr::rename(chromosome = Chrom,
			        position = Pos,
			        reference_allele = Ref,
			        total_depth = TOTAL_depth) %>%
		  dplyr::select(chromosome, position, reference_allele, total_depth, alternate_depth)
	write_tsv(pile_up, path = paste0("summary/", sample_name, "/ins_del/", bar_code, ".txt"), na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==3) {
	sample_name = opt$sample_name
	pile_up = readr::read_tsv(file = paste0("gbc/", sample_name, "/", bar_code, ".txt.gz"),
				  col_names = TRUE,
				  col_types = cols(.default = col_character())) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(A = A + a,
			        T = T + t,
			        G = G + g,
			        C = C + c) %>%
		  dplyr::rename(chromosome = Chrom,
			        position = Pos,
			        reference_allele = Ref,
			        total_depth = TOTAL_depth) %>%
		  dplyr::select(chromosome, position, reference_allele, total_depth, A, T, G, C) %>%
		  reshape2::melt(id.vars = c("chromosome", "position", "reference_allele", "total_depth"),
				 measure.vars = c("A", "T", "G", "C"),
				 variable.name = "alternate_allele",
				 value.name = "alternate_depth") %>%
		  dplyr::select(chromosome, position, reference_allele, alternate_allele, total_depth, alternate_depth) %>%
		  dplyr::filter(reference_allele != alternate_allele) %>%
		  dplyr::mutate(numeric_chromosome = gsub(pattern = "chr", replacement = "", x = chromosome, fixed = TRUE)) %>%
		  dplyr::mutate(numeric_chromosome = case_when(
			  numeric_chromosome == "X" ~ "23",
			  numeric_chromosome == "Y" ~ "24",
			  numeric_chromosome == "MT" ~ "25",
			  TRUE ~ numeric_chromosome
		  )) %>%
		  dplyr::mutate(numeric_chromosome = as.numeric(numeric_chromosome)) %>%
		  dplyr::arrange(numeric_chromosome, position, reference_allele, alternate_allele) %>%
		  dplyr::select(chromosome, position, reference_allele, alternate_allele, total_depth, alternate_depth)
	write_tsv(pile_up, path = paste0("summary/", sample_name, "/all_alt/", bar_code, ".txt"), na = "NA", append = FALSE, col_names = TRUE)

}
