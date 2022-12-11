#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("copynumber"))

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
		dplyr::select(Chromosome = chrom,
			      Position = maploc,
			      Log2_Ratio = cnlr)
	readr::write_tsv(x = cn_df, file = as.character(opt$file_out), col_names = TRUE, append = FALSE)

} else if (as.numeric(opt$option) == 2) {
	tumor_sample_names = unlist(strsplit(x = as.character(opt$tumor_sample_name), split = " ", fixed = TRUE))
	normal_sample_name = unlist(strsplit(x = as.character(opt$normal_sample_name), split = " ", fixed = TRUE))
	cn_df = list()
	for (i in 1:length(tumor_sample_names)) {
		data_ = readr::read_tsv(file = paste0("medicc2/", tumor_sample_names[i], "/", tumor_sample_names[i], ".txt"),
					col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert()
		colnames(data_) = c("Chromosome", "Position", paste0(tumor_sample_names[i], "_Log2_Ratio"))
		cn_df[[i]] = data_ %>%
			     reshape2::melt(id.vars = c("Chromosome", "Position"))
	}
	cn_df = do.call(bind_rows, cn_df) %>%
		reshape2::dcast(Chromosome + Position ~ variable, value.var = "value", fill = 0)
	readr::write_tsv(x = cn_df, file = as.character(opt$file_out), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option) == 3) {
	tumor_sample_names = unlist(strsplit(x = as.character(opt$tumor_sample_name), split = " ", fixed = TRUE))
	normal_sample_name = unlist(strsplit(x = as.character(opt$normal_sample_name), split = " ", fixed = TRUE))
	cn_df = readr::read_tsv(file = as.character(opt$file_in), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert() %>%
		as.data.frame()
	cn_smooth = copynumber::winsorize(data = cn_df, method = "mad", tau = 2.5, k = 25, verbose = FALSE)
	cn_segmented = copynumber::multipcf(data = cn_smooth, gamma = 40, normalize = FALSE, fast = FALSE, verbose = FALSE)
	
	total_copies = cn_segmented %>%
		       dplyr::select(c("chrom", "start.pos", "end.pos", contains("Log2_Ratio"))) %>%
		       dplyr::rename(start = start.pos, end = end.pos) %>%
		       reshape2::melt(id.vars = c("chrom", "start", "end")) %>%
		       dplyr::select(sample_id = variable,
				     chrom, start, end, nAB = value) %>%
		       dplyr::mutate(sample_id = gsub(pattern = "_Log2_Ratio", replacement = "", x = sample_id, fixed = TRUE)) %>%
		       dplyr::left_join(readr::read_tsv(file = "facets/summary/summary.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       				dplyr::select(sample_id = tumorName, purity, ploidy),
					by = "sample_id") %>%
		       readr::type_convert() %>%
		       dplyr::mutate(nAB = ((2^nAB)*((purity*ploidy) + (2*(1-purity))) - 2*(1-purity))/purity) %>%
		       dplyr::mutate(nAB = round(nAB)) %>%
		       dplyr::select(-purity, -ploidy)
	
	readr::write_tsv(x = total_copies, file = as.character(opt$file_out), col_names = TRUE, append = FALSE)
}
