#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("immunedeconv"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option('--option', type = 'character', default = NA, help = 'Immune deconv algorithm'),
               make_option('--input_file', type = 'character', default = NA, help = 'Expression input file'),
	       make_option('--output_file', type = 'character', default = NA, help = 'Immune cell output file'))
parser = OptionParser(usage = "%prog",  option_list=optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

set_cibersort_binary("~/share/usr/lib/resource_files/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("~/share/usr/lib/resource_files/CIBERSORT/LM22.txt")

if (as.numeric(opt$option)==1) {
	tpm_by_gene = readr::read_tsv(file = opt$input_file, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::arrange(gene_symbol)
	gene_expr = tpm_by_gene %>%
		    dplyr::select(-gene_symbol) %>%
		    as.matrix()
	rownames(gene_expr) = tpm_by_gene %>% .[["gene_symbol"]]
	quantiseq = immunedeconv::deconvolute(gene_expression = gene_expr, method = "quantiseq", scale_mrna = FALSE)
	readr::write_tsv(x = quantiseq, file = opt$output_file, col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==2) {
	tpm_by_gene = readr::read_tsv(file = opt$input_file, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::arrange(gene_symbol)
	gene_expr = tpm_by_gene %>%
		    dplyr::select(-gene_symbol) %>%
		    as.matrix()
	rownames(gene_expr) = tpm_by_gene %>% .[["gene_symbol"]]
	mcpcounter = immunedeconv::deconvolute(gene_expression = gene_expr, method = "mcp_counter", scale_mrna = FALSE)
	readr::write_tsv(x = mcpcounter, file = opt$output_file, col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==3) {
		tpm_by_gene = readr::read_tsv(file = opt$input_file, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::arrange(gene_symbol)
	gene_expr = tpm_by_gene %>%
		    dplyr::select(-gene_symbol) %>%
		    as.matrix()
	rownames(gene_expr) = tpm_by_gene %>% .[["gene_symbol"]]
	cibersort = immunedeconv::deconvolute(gene_expression = gene_expr, method = "cibersort_abs", scale_mrna = FALSE)
	readr::write_tsv(x = cibersort, file = opt$output_file, col_names = TRUE, append = FALSE)

}
