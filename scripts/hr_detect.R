#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "type of analysis"),
		  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option) == 1) {
	vcf = readr::read_tsv(file = "summary/tsv/all.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::mutate(TUMOR_NORMAL = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
	      dplyr::filter(TUMOR_NORMAL == as.character(opt$sample_name)) %>%
	      dplyr::filter(variantCaller == "mutect") %>%
	      dplyr::filter(Variant_Classification == "Missense_Mutation") %>%
	      dplyr::select(`#CHROM` = CHROM,
			    POS = POS,
			    ID = ID,
			    REF = REF,
			    ALT = ALT)
	cat("##fileformat=VCFv4.1\n", file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), append = FALSE)
	readr::write_tsv(x = vcf, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), col_names = TRUE, append = TRUE)

}
