#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1,
            error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "sample names"),
	       make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))

if (as.numeric(opt$option)==1) {
	
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		x1[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR-intervals.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
		  	  .[["X7"]]
	}
	x1 = do.call(cbind, x1)
	colnames(x1) = sample_names
	x1 = dplyr::as_tibble(x1) %>%
	     dplyr::mutate(LIBRARY = "STANDARD")

	for (i in 1:length(sample_names)) {
		x2[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-intervals.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		  	  readr::type_convert() %>%
			  .[["X7"]]
	}
	x2 = do.call(cbind, x2)
	colnames(x2) = sample_names
	x2 = dplyr::as_tibble(x2) %>%
	     dplyr::mutate(LIBRARY = "COLLAPSED")

	for (i in 1:length(sample_names)) {
		x3[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-intervals.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  .[["X7"]]
	}
	x3 = do.call(cbind, x3)
	colnames(x3) = sample_names
	x3 = dplyr::as_tibble(x3) %>%
	     dplyr::mutate(LIBRARY = "SIMPLEX")

	for (i in 1:length(sample_names)) {
		x4[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-intervals.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  .[["X7"]]
	}
	x4 = do.call(cbind, x4)
	colnames(x4) = sample_names
	x4 = dplyr::as_tibble(x4) %>%
	     dplyr::mutate(LIBRARY = "DUPLEX")

	x = x1 %>%
	bind_rows(x2) %>%
	bind_rows(x3) %>%
	bind_rows(x4)
	write_tsv(x, path="summary/intervals_summary.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==2) {
	
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		x1[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x1 = do.call(rbind, x1)
	x1 = dplyr::as_tibble(x1) %>%
	     dplyr::mutate(LIBRARY = "STANDARD") %>%
	     dplyr::rename(Chromosome = X1,
			   Position = X2,
			   Reference_Allele = X3,
			   Total_Count = X4,
			   A = X5,
			   C = X6,
			   G = X7,
			   T = X8,
			   INS = X9,
			   DEL = X10,
			   Soft_Clip_Start = X11,
			   Soft_Clip_End = X12,
			   Hard_Clip_Start = X13,
			   Hard_Clip_End = X14) %>%
	     dplyr::mutate(N = Total_Count - (A + C + G + T))
	
	for (i in 1:length(sample_names)) {
		x2[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		  	  readr::type_convert() %>%
			  dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x2 = do.call(rbind, x2)
	x2 = dplyr::as_tibble(x2) %>%
	     dplyr::mutate(LIBRARY = "COLLAPSED") %>%
	     dplyr::rename(Chromosome = X1,
			   Position = X2,
			   Reference_Allele = X3,
			   Total_Count = X4,
			   A = X5,
			   C = X6,
			   G = X7,
			   T = X8,
			   INS = X9,
			   DEL = X10,
			   Soft_Clip_Start = X11,
			   Soft_Clip_End = X12,
			   Hard_Clip_Start = X13,
			   Hard_Clip_End = X14) %>%
	     dplyr::mutate(N = Total_Count - (A + C + G + T))
	
	for (i in 1:length(sample_names)) {
		x3[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x3 = do.call(rbind, x3)
	x3 = dplyr::as_tibble(x3) %>%
	     dplyr::mutate(LIBRARY = "SIMPLEX") %>%
	     dplyr::rename(Chromosome = X1,
			   Position = X2,
			   Reference_Allele = X3,
			   Total_Count = X4,
			   A = X5,
			   C = X6,
			   G = X7,
			   T = X8,
			   INS = X9,
			   DEL = X10,
			   Soft_Clip_Start = X11,
			   Soft_Clip_End = X12,
			   Hard_Clip_Start = X13,
			   Hard_Clip_End = X14) %>%
	     dplyr::mutate(N = Total_Count - (A + C + G + T))
	
	for (i in 1:length(sample_names)) {
		x4[[i]] = readr::read_tsv(file=paste0("waltz/", sample_names[i], "_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  dplyr::mutate(SAMPLE_NAME = sample_names[i])
	}
	x4 = do.call(rbind, x4)
	x4 = dplyr::as_tibble(x4) %>%
	     dplyr::mutate(LIBRARY = "DUPLEX") %>%
	     dplyr::rename(Chromosome = X1,
			   Position = X2,
			   Reference_Allele = X3,
			   Total_Count = X4,
			   A = X5,
			   C = X6,
			   G = X7,
			   T = X8,
			   INS = X9,
			   DEL = X10,
			   Soft_Clip_Start = X11,
			   Soft_Clip_End = X12,
			   Hard_Clip_Start = X13,
			   Hard_Clip_End = X14) %>%
	     dplyr::mutate(N = Total_Count - (A + C + G + T))

	write_tsv(x1, path="summary/pileup_summary_standard.txt", na = "NA", append = FALSE, col_names = TRUE)
	write_tsv(x2, path="summary/pileup_summary_collapsed.txt", na = "NA", append = FALSE, col_names = TRUE)
	write_tsv(x3, path="summary/pileup_summary_simplex.txt", na = "NA", append = FALSE, col_names = TRUE)
	write_tsv(x4, path="summary/pileup_summary_duplex.txt", na = "NA", append = FALSE, col_names = TRUE)
}
