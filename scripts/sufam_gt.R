#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("fuzzyjoin"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_set", default = NA, type = 'character', help = "sample set"),
	       make_option("--tumor_sample", default = NA, type = 'character', help = "tumor sample"),
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
	       dplyr::mutate(`#CHROM` = as.character(`#CHROM`)) %>%
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
	tumor_sample = unlist(strsplit(x = as.character(opt$tumor_sample), split = " ", fixed=TRUE))
	normal_sample = unlist(strsplit(x = as.character(opt$normal_sample), split = " ", fixed=TRUE))
	maf = readr::read_tsv(file = opt$input_file, comment = "#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::mutate(chrom = Chromosome,
			    loc.start = Start_Position,
			    loc.end = End_Position) %>%
	      dplyr::mutate(chrom = as.character(chrom))
	facets = readr::read_tsv(file = paste0("facets/cncf/", tumor_sample, "_", normal_sample, ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 dplyr::mutate(chrom = case_when(
			 chrom == "23" ~ "X",
			 TRUE ~ chrom
		 )) %>%
		 readr::type_convert() %>%
		 dplyr::mutate(qt = tcn.em,
			       q2 = tcn.em - lcn.em) %>%
		 dplyr::select(chrom, loc.start, loc.end, qt, q2)
	maf = maf %>%
	      fuzzyjoin::genome_left_join(facets, by = c("chrom", "loc.start", "loc.end")) %>%
	      dplyr::select(-chrom.x, -loc.start.x, -loc.end.x, -chrom.y, -loc.start.y, -loc.end.y)
			
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==99) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	normal_sample = unlist(strsplit(x = as.character(opt$normal_sample), split = " ", fixed=TRUE))
	sample_set = setdiff(sample_set, normal_sample)
	maf = list()
	for (i in 1:length(sample_set)) {
		sufam = readr::read_tsv(file = paste0("sufam/", sample_set[i], ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(CHROM = chrom,
				      POS = pos,
				      REF = val_ref,
				      ALT = val_alt,
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
	      dplyr::filter(t_alt_count > 1)
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==4) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	maf = list()
	for (i in 1:length(sample_set)) {
		maf[[i]] = readr::read_tsv(file = paste0("sufam/", sample_set[i], ".maf"), comment = "#", col_names = TRUE, col_types = cols(.default = col_character()))
	}
	maf = do.call(bind_rows, maf) %>%
	readr::type_convert()
	smry = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       dplyr::mutate(HOTSPOT = case_when(
		       is.na(HOTSPOT) ~ FALSE,
		       HOTSPOT == "True" ~ TRUE,
		       HOTSPOT == "False" ~ FALSE,
		       HOTSPOT == "TRUE" ~ TRUE,
		       HOTSPOT == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(HOTSPOT_INTERNAL = case_when(
		       is.na(HOTSPOT_INTERNAL) ~ FALSE,
		       HOTSPOT_INTERNAL == "True" ~ TRUE,
		       HOTSPOT_INTERNAL == "False" ~ FALSE,
		       HOTSPOT_INTERNAL == "TRUE" ~ TRUE,
		       HOTSPOT_INTERNAL == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(cmo_hotspot = case_when(
		       is.na(cmo_hotspot) ~ FALSE,
		       cmo_hotspot == "True" ~ TRUE,
		       cmo_hotspot == "False" ~ FALSE,
		       cmo_hotspot == "TRUE" ~ TRUE,
		       cmo_hotspot == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(is_hotspot = HOTSPOT | HOTSPOT_INTERNAL | cmo_hotspot) %>%
	       dplyr::mutate(facetsLOHCall = case_when(
		       is.na(facetsLOHCall) ~ FALSE,
		       facetsLOHCall == "True" ~ TRUE,
		       facetsLOHCall == "False" ~ FALSE,
		       facetsLOHCall == "TRUE" ~ TRUE,
		       facetsLOHCall == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(is_loh = facetsLOHCall) %>%
	       readr::type_convert()
	maf = maf %>%
	      dplyr::left_join(smry %>%
			       dplyr::group_by(CHROM, POS, REF, ALT) %>%
	       		       dplyr::summarize(is_hotspot = unique(is_hotspot)) %>%
			       dplyr::ungroup(),
			       by = c("CHROM", "POS", "REF", "ALT"))
	maf = maf %>%
	      dplyr::left_join(smry %>%
			       dplyr::select(CHROM, POS, REF, ALT, Tumor_Sample_Barcode = TUMOR_SAMPLE, is_loh),
			       by = c("CHROM", "POS", "REF", "ALT", "Tumor_Sample_Barcode"))
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$option)==5) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " ", fixed=TRUE))
	maf = list()
	for (i in 1:length(sample_set)) {
		maf[[i]] = readr::read_tsv(file = paste0("sufam/", sample_set[i], "_ft.maf"), comment = "#", col_names = TRUE, col_types = cols(.default = col_character()))
	}
	maf = do.call(bind_rows, maf) %>%
	readr::type_convert()
	smry = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       dplyr::mutate(HOTSPOT = case_when(
		       is.na(HOTSPOT) ~ FALSE,
		       HOTSPOT == "True" ~ TRUE,
		       HOTSPOT == "False" ~ FALSE,
		       HOTSPOT == "TRUE" ~ TRUE,
		       HOTSPOT == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(HOTSPOT_INTERNAL = case_when(
		       is.na(HOTSPOT_INTERNAL) ~ FALSE,
		       HOTSPOT_INTERNAL == "True" ~ TRUE,
		       HOTSPOT_INTERNAL == "False" ~ FALSE,
		       HOTSPOT_INTERNAL == "TRUE" ~ TRUE,
		       HOTSPOT_INTERNAL == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(cmo_hotspot = case_when(
		       is.na(cmo_hotspot) ~ FALSE,
		       cmo_hotspot == "True" ~ TRUE,
		       cmo_hotspot == "False" ~ FALSE,
		       cmo_hotspot == "TRUE" ~ TRUE,
		       cmo_hotspot == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(is_hotspot = HOTSPOT | HOTSPOT_INTERNAL | cmo_hotspot) %>%
	       dplyr::mutate(facetsLOHCall = case_when(
		       is.na(facetsLOHCall) ~ FALSE,
		       facetsLOHCall == "True" ~ TRUE,
		       facetsLOHCall == "False" ~ FALSE,
		       facetsLOHCall == "TRUE" ~ TRUE,
		       facetsLOHCall == "FALSE" ~ FALSE
	       )) %>%
	       dplyr::mutate(is_loh = facetsLOHCall) %>%
	       readr::type_convert()
	maf = maf %>%
	      dplyr::left_join(smry %>%
			       dplyr::group_by(CHROM, POS, REF, ALT) %>%
	       		       dplyr::summarize(is_hotspot = unique(is_hotspot)) %>%
			       dplyr::ungroup(),
			       by = c("CHROM", "POS", "REF", "ALT"))
	maf = maf %>%
	      dplyr::left_join(smry %>%
			       dplyr::select(CHROM, POS, REF, ALT, Tumor_Sample_Barcode = TUMOR_SAMPLE, is_loh),
			       by = c("CHROM", "POS", "REF", "ALT", "Tumor_Sample_Barcode"))
	write_tsv(x = maf, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)

}

