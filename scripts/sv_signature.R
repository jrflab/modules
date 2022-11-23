#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("signature.tools.lib"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_name = as.character(opt$sample_name)
	sv_bedpe = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   readr::type_convert() %>%
		   dplyr::mutate(sample = sample_name) %>%
		   dplyr::select(-svclass)
	res_list = bedpeToRearrCatalogue(sv_bedpe %>% data.frame())
	catalogues_mutations = data.frame(row.names = rownames(res_list$rearr_catalogue), stringsAsFactors = FALSE)
	bedpecolumns = c("chrom1", "start1", "end1", "chrom2", "start2", "end2" , "sample","svclass","id", "is.clustered", "length")
	catalogues_mutations = cbind(catalogues_mutations,res_list$rearr_catalogue)
	mtype_mutations = signature.tools.lib:::getTypeOfMutationsFromChannels(catalogues_mutations)
	exposureFilterType = "fixedThreshold"
	threshold_percent = 5
	optimisation_method = "KLD"
	useBootstrap = FALSE
	nboot = 1000
	threshold_p.value = 0.05
	nparallel = 4
	randomSeed = 1
	fit = Fit(catalogues = catalogues_mutations,
		  signatures = signature.tools.lib:::RefSigv1_rearr,
		  exposureFilterType = exposureFilterType,
		  threshold_percent = threshold_percent,
		  method = optimisation_method,
		  useBootstrap = useBootstrap,
		  nboot = nboot,
		  threshold_p.value = threshold_p.value,
		  nparallel = nparallel,
		  randomSeed = randomSeed,
		  verbose = TRUE)
	x = dplyr::tibble(feature_name = rownames(fit$catalogues),
			  feature_count = as.vector(fit$catalogues[,1])) %>%
	    dplyr::mutate(sample_name = sample_name)
	readr::write_tsv(x = x, file = paste0(opt$output_file, "_features.txt"), col_names = TRUE, append = FALSE)
	
	x = dplyr::tibble(signature_name = colnames(fit$exposures),
			  signature_exposure = as.vector(fit$exposures[1,])/sum(as.vector(fit$exposures[1,])) * 100) %>%
	    dplyr::mutate(sample_name = sample_name)
	readr::write_tsv(x = x, file = paste0(opt$output_file, "_exposures.txt"), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==2) {
	sample_name = as.character(opt$sample_name)
	bedpe_org = readr::read_tsv(file = paste0("sv_signature/", sample_name, "/", sample_name, ".merged.bedpe"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    dplyr::filter(chrom1 != "Y") %>%
	    	    dplyr::filter(chrom2 != "Y") %>%
		    readr::type_convert()
	bedpe_cli = readr::read_tsv(file = paste0("sv_signature/", sample_name, "/", sample_name, ".merged.sv_clusters_and_footprints.tsv"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert() %>%
		    dplyr::select(chrom1 = X1,
				  start1 = X2,
				  end1 = X3,
				  chrom2 = X4,
				  start2 = X5,
				  end2 = X6,
				  n_svs = X12,
				  p_value = X17) %>%
		    dplyr::mutate(p_value = as.numeric(p_value)) %>%
		    dplyr::mutate(p_value = case_when(
			  is.na(p_value) ~ 1,
			  TRUE ~ p_value))
	bedpe_org = bedpe_org %>%
		    dplyr::left_join(bedpe_cli, by = c("chrom1", "start1", "end1", "chrom2", "start2", "end2")) %>%
	    	    dplyr::mutate(is_clustered = case_when(
		    	p_value<.05 & n_svs>=15 ~ "c1",
		    	TRUE ~ "non_clustered"
	    	    )) %>%
		    dplyr::mutate(is_clustered = case_when(
		    	p_value<.05 & n_svs>=50 ~ "c2",
		    	TRUE ~ is_clustered
		    )) %>%
		    dplyr::mutate(svclass = case_when(
		    	svclass == "TRA" & is_clustered == "c1" ~ "c1TRA",
		    	svclass == "TRA" & is_clustered == "c2" ~ "c2TRA",
		    	svclass == "INV" & (is_clustered == "c1" | is_clustered == "c2") ~ "cINV",
		    	TRUE ~ svclass
		    )) %>%
		    dplyr::select(chrom1, start1, end1, chrom2, start2, end2, sv_id, pe_support, strand1, strand2, svclass)
	write_tsv(x = bedpe_org, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==3) {
	sample_name = as.character(opt$sample_name)
	catalogues = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		     readr::type_convert()
	catalogues_mutations = data.frame(catalogues %>% dplyr::select(manual_sv_type))
	colnames(catalogues_mutations) = sample_name
	rownames(catalogues_mutations) = catalogues %>% .[["...1"]]
	
	signatures = readr::read_tsv(file = "~/share/lib/resource_files/viola/NMF/signature_matrix.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		     readr::type_convert()
	signatures_mutations = data.frame(signatures %>% dplyr::select(-`SV Type`))
	colnames(signatures_mutations) = colnames(signatures)[-1]
	rownames(signatures_mutations) = signatures %>% .[["SV Type"]]
	exposureFilterType = "fixedThreshold"
	threshold_percent = 5
	optimisation_method = "KLD"
	useBootstrap = FALSE
	nboot = 1000
	threshold_p.value = 0.05
	nparallel = 4
	randomSeed = 1
	fit = Fit(catalogues = catalogues_mutations,
		  signatures = signatures_mutations,
		  exposureFilterType = exposureFilterType,
		  threshold_percent = threshold_percent,
		  method = optimisation_method,
		  useBootstrap = useBootstrap,
		  nboot = nboot,
		  threshold_p.value = threshold_p.value,
		  nparallel = nparallel,
		  randomSeed = randomSeed,
		  verbose = TRUE)
	x = dplyr::tibble(feature_name = rownames(fit$catalogues),
			  feature_count = as.vector(fit$catalogues[,1])) %>%
	    dplyr::mutate(sample_name = sample_name)
	readr::write_tsv(x = x, file = paste0(opt$output_file, "_features.txt"), col_names = TRUE, append = FALSE)
	
	x = dplyr::tibble(signature_name = colnames(fit$exposures),
			  signature_exposure = as.vector(fit$exposures[1,])/sum(as.vector(fit$exposures[1,])) * 100) %>%
	    dplyr::mutate(sample_name = sample_name)
	readr::write_tsv(x = x, file = paste0(opt$output_file, "_exposures.txt"), col_names = TRUE, append = FALSE)

}
