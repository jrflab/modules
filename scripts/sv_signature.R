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
	       make_option("--input_file", default = NA, type = 'character', help = "input file")
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
	readr::write_tsv(x = fit$catalogues, file = paste0(opt$output_file, "_features.txt"), col_names = TRUE, append = FALSE)
	readr::write_tsv(x = fit$catalogues, file = paste0(opt$output_file, "_exposures.txt"), col_names = TRUE, append = FALSE)
	
}

#if (as.numeric(opt$option)==1) {
#	sample_names = as.character(opt$sample_names)
#	bedpe_org = readr::read_tsv(file = paste0("sv_signature/", sample_names, "/", sample_names, ".merged.bedpe"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
#		    readr::type_convert()
#	bedpe_cli = readr::read_tsv(file = paste0("sv_signature/", sample_names, "/", sample_names, ".sv_clusters_and_footprints.tsv"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
#		    readr::type_convert() %>%
#		    dplyr::select(chrom1 = X1,
#				  start1 = X2,
#				  end1 = X3,
#				  chrom2 = X4,
#				  start2 = X5,
#				  end2 = X6,
#				  sv_id = X7,
#				  pe_support = X8,
#				  strand1 = X9,
#				  strand2 = X10,
#				  n_svs = X12,
#				  p_value = X17) %>%
#		    dplyr::mutate(p_value = case_when(
#			    is.na(p_value) ~ 1,
#			    TRUE ~ p_value
#		    )) %>%
#		    dplyr::left_join(bedpe_org %>%
#				     dplyr::select(chrom1, start1, end1, chrom2, start2, end2, svclass),
#				     by = c("chrom1", "start1", "end1", "chrom2", "start2", "end2")) %>%
#		    dplyr::filter(!is.na(svclass)) %>%
#		    dplyr::mutate(is_clustered = case_when(
#			    p_value < as.numeric(opt$p_value) & n_svs > 2*as.numeric(opt$n_sv) ~ "Cplx2",
#			    p_value < as.numeric(opt$p_value) & n_svs > as.numeric(opt$n_sv) ~ "Cplx1",
#			    TRUE ~ ""
#		    )) %>%
#		    dplyr::mutate(svclass = case_when(
#			    svclass == "TRA" & is_clustered != "" ~ paste0(is_clustered, svclass),
#			    TRUE ~ svclass
#		    )) %>%
#		    dplyr::select(chrom1, start1, end1, chrom2, start2, end2, sv_id, pe_support, strand1, strand2, svclass)
#	
#	write_tsv(x = bedpe_cli, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
#	
#} else if (as.numeric(opt$option)==2) {
#	sample_names = unlist(strsplit(x = as.character(opt$sample_names), split = " ", fixed=TRUE))
#	feature_counts = list()
#	for (i in 1:length(sample_names)) {
#		feature_counts[[i]] = readr::read_tsv(file = paste0("sv_signature/", sample_names[i], "/", sample_names[i], ".merged.txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
#				      readr::type_convert() %>%
#				      dplyr::rename(sv_class = X1,
#					            sv_count = manual_sv_type) %>%
#				      dplyr::mutate(sample_name = sample_names[i])
#	}
#	feature_counts = do.call(bind_rows, feature_counts)
#	write_tsv(x = feature_counts, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
#	
#}
