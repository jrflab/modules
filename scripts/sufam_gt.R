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
	       make_option("--normal_sample", default = NA, type = 'character', help = "normal sample"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	normal_sample = unlist(strsplit(x = as.character(opt$normal_sample), split = " ", fixed=TRUE))
	sample_set = setdiff(sample_set, normal_sample)
	smry = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::filter(TUMOR_SAMPLE %in% sample_set) %>%
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

} else if (as.numeric(opt$option)==2) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	normal_sample = unlist(strsplit(x = as.character(opt$normal_sample), split = " ", fixed=TRUE))
	sample_set = setdiff(sample_set, normal_sample)
	maf = list()
	for (i in 1:length(sample_set)) {
		sufam = readr::read_tsv(file = paste0("sufam/", sample_set[i], ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(CHROM = chrom,
				      POS = pos,
				      REF = ref,
				      t_depth = cov,
				      t_alt_count = val_al_count) %>%
		 	dplyr::mutate(t_ref_count = t_depth - t_alt_count)
			
		maf[[i]] = readr::read_tsv(file = paste0("sufam/", sample_set[i], ".maf"), comment = "#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      	   readr::type_convert() %>%
		      	   dplyr::select(-t_depth, -t_alt_count, -t_ref_count) %>%
		      	   dplyr::bind_cols(sufam)
	}
	maf = do.call(bind_rows, maf)
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==3) {
	maf = readr::read_tsv(file = as.character(opt$input_file), comment = "#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::filter(t_alt_count > 0) %>%
	      dplyr::filter(t_ref_count > 0)
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==4) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	maf = list()
	for (i in 1:length(sample_set)) {
		maf[[i]] = readr::read_tsv(file = paste0("sufam/", sample_set[i], ".maf"), comment = "#", col_names = TRUE, col_types = cols(.default = col_character()))
	}
	maf = do.call(bind_rows, maf)
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==5) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	maf = list()
	for (i in 1:length(sample_set)) {
		maf[[i]] = readr::read_tsv(file = paste0("sufam/", sample_set[i], "_ft.maf"), comment = "#", col_names = TRUE, col_types = cols(.default = col_character()))
	}
	maf = do.call(bind_rows, maf)
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

}

