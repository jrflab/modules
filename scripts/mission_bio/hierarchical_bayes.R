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
	nb = 5E3
	ni = 1E3
	data = list(yp = pile_up %>% .[["alternate_depth"]],
		    dp = pile_up %>% dplyr::mutate(total_depth = log(total_depth)) %>% .[["total_depth"]],
	    	    n3 = pile_up %>% .[["levels"]],
	    	    N = nrow(pile_up),
	    	    L = max(pile_up %>% .[["levels"]]))

	params = c('lambda', 'alpha', 'a', 'b', 'tau', 'tau.b', 'psi')
	post_mcmc = ZIPR(data = data,
			 params = params,
			 target = paste0("hbm/", sample_name, "/mcmc/", bar_code, ".jags"),
			 nc = nc, nd = nd, nb = nb, ni = ni)
	save(list = ls(all = TRUE), file = paste0("hbm/", sample_name, "/mcmc/", bar_code, ".RData"))

} else if (as.numeric(opt$option)==2) {
	sample_name = opt$sample_name
	bar_code = opt$bar_code
	
	load(file = paste0("hbm/", sample_name, "/mcmc/", bar_code, ".RData"))
	
	post_est = apply(do.call(rbind, post_mcmc), 2, median)
	index_lambda = grep("lambda[", names(post_est), fixed = TRUE)
	index_alpha = which(names(post_est) == "alpha")
	index_a = grep("^a\\[", names(post_est), perl = TRUE)
	index_b = grep("^b\\[", names(post_est), perl = TRUE)
	index_tau = which(names(post_est) == "tau")
	index_tau.b = which(names(post_est) == "tau.b")
	index_psi = which(names(post_est) == "psi")
	
	lambda_p = post_est[index_lambda]
	alpha = post_est[index_alpha]
	a = post_est[index_a]
	b = post_est[index_b]
	tau = post_est[index_tau]
	tau_b = post_est[index_tau.b]
	psi = post_est[index_psi]
	
	nucleotide_context = pile_up %>%
			     dplyr::group_by(context_3) %>%
			     dplyr::summarize(level = unique(levels)) %>%
			     dplyr::ungroup() %>%
			     dplyr::arrange(level) %>%
			     dplyr::mutate(a = as.vector(a)) %>%
			     dplyr::select(-level) %>%
			     dplyr::rename(context = context_3)
	
	a = nucleotide_context %>% .[["a"]]
	names(a) = nucleotide_context %>% .[["context"]]
	
	data = pile_up %>%
	       dplyr::select(chromosome, position, context_3) %>%
	       dplyr::mutate(y_p = data$yp, d_p = data$dp) %>%
	       dplyr::mutate(d_p = round(exp(d_p))) %>%
	       dplyr::mutate(lambda_p = lambda_p, b_p = b) %>%
	       dplyr::mutate(alpha = alpha, tau = tau, tau_b = tau_b, psi = psi) %>%
	       dplyr::mutate(a = a[context_3])
	
	save(data, file = paste0("hbm/", sample_name, "/params/", bar_code, ".RData"))
	
} else if (as.numeric(opt$option)==3) {
	sample_name = as.character(opt$sample_name)
	bar_code = unlist(strsplit(x = as.character(opt$bar_code), split = " ", fixed = TRUE))

	lambda_p = d_p = list()
	for (i in 1:length(bar_code)) {
		load(file = paste0("hbm/", sample_name, "/params/", bar_code[i], ".RData"))
		#data = data %>%
		#       dplyr::mutate(lambdap_dp = lambda_p + tau_b*b_p) %>%
		#       dplyr::mutate(lambdap_dp = case_when(
		#	       d_p == 0 ~ 100*lambdap_dp,
		#	       d_p > 0 ~ 100*(lambdap_dp+2)/(d_p+4)
		#       )) %>%
		#       dplyr::select(lambdap_dp)
		data = data %>%
		       dplyr::mutate(lambda_p = lambda_p + tau_b*b_p)
		lambda_p[[i]] = data$lambda_p
		d_p[[i]] = data$d_p
	}
	lambda_p = do.call(bind_cols, lambda_p) %>%
		   dplyr::as_tibble()
	d_p = do.call(bind_cols, d_p) %>%
	      dplyr::as_tibble()
	colnames(lambda_p) = colnames(d_p) = bar_code
	
	load(file = paste0("hbm/", sample_name, "/params/", bar_code[i], ".RData"))
	lambda_p = dplyr::bind_cols(data %>% dplyr::select(chromosome, position), lambda_p)
	d_p = dplyr::bind_cols(data %>% dplyr::select(chromosome, position), d_p)
	save(lambda_p, d_p, file = paste0("hbm/", sample_name, "/__LambdaP_dP__.RData"))

}
