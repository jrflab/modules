#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("MSKscmbio"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
	       make_option("--snp_file", default = NA, type = 'character', help = "germline SNP"),
	       make_option("--context_file", default = NA, type = 'character', help = "nucleotide context"),
               make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
	       make_option("--bar_code", default = NA, type = 'character', help = "sample name"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_name = opt$sample_name
	bar_code = opt$bar_code
	
	germline_snp = readr::read_tsv(file = opt$snp_file, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		       readr::type_convert() %>%
		       dplyr::mutate(is_snp = TRUE)
	
	nucleotide_context = readr::read_tsv(file = opt$context_file, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			     dplyr::group_by(chromosome, position, ref) %>%
	      		     dplyr::summarize(left_1nt = unique(left_1nt),
					      right_1nt = unique(right_1nt),
					      left_2nt = unique(left_2nt),
					      right_2nt = unique(right_2nt)) %>%
			     dplyr::ungroup() %>%
			     dplyr::mutate(context_3 = paste0(left_1nt, ref, right_1nt)) %>%
			     dplyr::mutate(context_5 = paste0(left_2nt, ref, right_2nt)) %>%
			     readr::type_convert()
	
	pile_up = readr::read_tsv(file = paste0("summary/", sample_name, "/sum_alt/", bar_code, ".txt.gz"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		  readr::type_convert() %>%
		  dplyr::left_join(germline_snp, by = c("chromosome", "position", "reference_allele")) %>%
		  dplyr::mutate(is_snp = case_when(
			is.na(is_snp) ~ FALSE,
			TRUE ~TRUE
		  )) %>%
		  dplyr::filter(!is_snp) %>%
		  dplyr::left_join(nucleotide_context %>% dplyr::select(chromosome, position, reference_allele = ref, context_3),
				   by = c("chromosome", "position", "reference_allele")) %>%
		  dplyr::select(-alternate_allele, -is_snp) %>%
		  dplyr::left_join(dplyr::tibble(context_3 = unique(nucleotide_context %>% .[["context_3"]])) %>%
		      		   dplyr::mutate(levels  = 1:nrow(.)), by = "context_3")
	
	#################################################
	## zero-inflated poisson regression
	#################################################
	nc = 2
	nd = 1E3
	nb = 1E4
	ni = 1E4
	data = list(yp = pile_up %>% .[["alternate_depth"]],	
		    dp = pile_up %>% .[["total_depth"]],
	    	    n3 = pile_up %>% .[["levels"]],
	    	    N = nrow(pile_up),
	    	    L = max(pile_up %>% .[["levels"]]))

	params = c('lambda', 'b', 'tau.b')
	post_mcmc = ZIPR(data = data,
			 params = params,
			 target = paste0("hierarchical_bayes/", sample_name, "/", bar_code, ".jags"),
			 nc = nc, nd = nd, nb = nb, ni = ni)
	save(list = ls(all = TRUE), paste0("hierarchical_bayes/", sample_name, "/", bar_code, ".RData"))
		   
}
