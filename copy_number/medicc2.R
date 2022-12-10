#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("reshape2"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "Which option?"),
		  make_option("--tumor_sample_name", default = NA, type = 'character', help = "Tumor sample name"),
		  make_option("--normal_sample_name", default = NA, type = 'character', help = "Normal sample name"),
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
			      B_Allele_F = 1 - vafT)
	readr::write_tsv(x = cn_df, file = as.character(opt$file_out), col_names = TRUE, append = FALSE)

} else if (as.numeric(opt$option) == 2) {
	tumor_sample_names = unlist(strsplit(x = as.character(opt$tumor_sample_name), split = " ", fixed = TRUE))
	normal_sample_name = unlist(strsplit(x = as.character(opt$normal_sample_name), split = " ", fixed = TRUE))
	cn_df = list()
	for (i in 1:length(tumor_sample_names)) {
		data_ = readr::read_tsv(file = paste0("medicc2/", tumor_sample_names[i], "/", tumor_sample_names[i], ".txt"),
					col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert()
		colnames(data_) = c("Chromosome", "Position", paste0(tumor_sample_names[i], "_Log2_Ratio"), paste0(tumor_sample_names[i], "_B_Allele_F"))
		cn_df[[i]] = data_ %>%
			     reshape2::melt(id.vars = c("Chromosome", "Position"))
	}
	cn_df = do.call(bind_rows, cn_df) %>%
		reshape2::dcast(Chromosome + Position ~ variable, value.var = "value", fill = 0)
	readr::write_tsv(x = cn_df, file = as.character(opt$file_out), col_names = TRUE, append = FALSE)
	
}
