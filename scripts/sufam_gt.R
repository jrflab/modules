#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_set", default = NA, type = 'character', help = "sample set"),
	       make_option("--normal_samples", default = NA, type = 'character', help = "normal samples"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_set), split="_", fixed=TRUE))
	normal_sample = intersect(sample_names, unlist(strsplit(x = as.character(opt$normal_samples), split=" ", fixed=TRUE)))
	sample_names = setdiff(sample_names, normal_sample)
	smry = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::filter(TUMOR_SAMPLE %in% sample_names) %>%
	       dplyr::filter(NORMAL_SAMPLE == normal_sample) %>%
	       dplyr::mutate(UUID = paste0(CHROM, ":", POS, "_", REF, ">", ALT)) %>%
	       dplyr::filter(!duplicated(UUID)) %>%
	       dplyr::mutate(`#CHROM` = CHROM,
			     POS = POS,
			     ID = ".",
			     REF = REF,
			     ALT = ALT,
			     QUAL = 100,
			     FILTER = "PASS",
			     INFO = ".") %>%
	       dplyr::select(`#CHROM`, POS, ID, REF, ALT, QUAL, INFO) %>%
	       dplyr::mutate(chr_n = case_when(
		       `#CHROM` == "X" ~ "23",
		       `#CHROM` == "Y" ~ "24",
		       TRUE ~ `#CHROM`
	       )) %>%
	       readr::type_convert() %>%
	       dplyr::arrange(chr_n) %>%
	       dplyr::select(-chr_n)
	
	cat("##fileformat=VCFv4.2\n", file = as.character(opt$output_file), append=FALSE)
	readr::write_tsv(x = smry, path = as.character(opt$output_file), append = TRUE, col_names = TRUE)
}

#all_vars = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#CHROM = all_vars[,"CHROM"]
#POS = all_vars[,"POS"]
#ID = all_vars[,"ID"]
#REF = all_vars[,"REF"]
#ALT = all_vars[,"ALT"]
#QUAL = FILTER = rep(".", nrow(all_vars))
#INFO = paste0(all_vars[,"SYMBOL"], all_vars[,"HGVSp_Short"])
#vcf = data.frame(CHROM, POS, ID, REF, ALT, QUAL, INFO)
#
#cat("#", file="sufam/pdx.vcf", append=FALSE)
#write.table(vcf, file="sufam/pdx.vcf", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
