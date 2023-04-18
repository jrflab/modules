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
	      dplyr::filter(TUMOR_DP>=10 & NORMAL_DP>=10) %>%
	      dplyr::mutate(CHROM = as.character(CHROM)) %>%
	      dplyr::mutate(CHROM = ifelse(CHROM == "23", "X", CHROM)) %>%
	      dplyr::select(`#CHROM` = CHROM,
			    POS = POS,
			    ID = ID,
			    REF = REF,
			    ALT = ALT)
	cat("##fileformat=VCFv4.1\n", file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), append = FALSE)
	readr::write_tsv(x = vcf, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), col_names = TRUE, append = TRUE)

} else if (as.numeric(opt$option) == 2) {
	vcf = readr::read_tsv(file = "summary/tsv/all.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::mutate(TUMOR_NORMAL = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
	      dplyr::filter(TUMOR_NORMAL == as.character(opt$sample_name)) %>%
	      dplyr::filter(grepl("varscan", variantCaller, fixed = TRUE)) %>%
	      dplyr::filter(grepl("strelka", variantCaller, fixed = TRUE)) %>%
	      dplyr::filter(TUMOR_DP>=10 & NORMAL_DP>=10) %>%
	      dplyr::mutate(CHROM = as.character(CHROM)) %>%
	      dplyr::mutate(CHROM = ifelse(CHROM == "23", "X", CHROM)) %>%
	      dplyr::select(`#CHROM` = CHROM,
			    POS = POS,
			    ID = ID,
			    REF = REF,
			    ALT = ALT)
	cat("##fileformat=VCFv4.1\n", file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".indel.vcf"), append = FALSE)
	readr::write_tsv(x = vcf, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".indel.vcf"), col_names = TRUE, append = TRUE)

} else if (as.numeric(opt$option) == 3) {
	cn = readr::read_tsv(file = paste0("facets/cncf/", as.character(opt$sample_name), ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
     	     dplyr::mutate(chrom = as.character(chrom)) %>%
	     dplyr::mutate(chrom = ifelse(chrom == "23", "X", chrom)) %>%
	     dplyr::mutate(seg_no = seg,
			   Chromosome = chrom,
			   chromStart = loc.start,
			   chromEnd = loc.end,
			   total.copy.number.inNormal = 2,
			   minor.copy.number.inNormal = 1,
			   total.copy.number.inTumour = tcn.em,
			   minor.copy.number.inTumour = lcn.em) %>%
	     dplyr::select(seg_no,
			   Chromosome,
			   chromStart,
			   chromEnd,
			   total.copy.number.inNormal,
			   minor.copy.number.inNormal,
			   total.copy.number.inTumour,
			   minor.copy.number.inTumour)
	     
	readr::write_tsv(x = cn, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".cn.txt"), col_names = TRUE, append = TRUE)
}
