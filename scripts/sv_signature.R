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
	sample_name = unlist(strsplit(x = as.character(opt$sample_name), split = " ", fixed = TRUE))
	signature_df = list()
	for (i in 1:length(sample_name)) {
		signature_df[[i]] = readr::read_tsv(file = paste0("sv_signature/", sample_name[i], "/", sample_name[i], ".merged_exposures.txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
				    readr::type_convert()
	}
	signature_df = do.call(bind_rows, signature_df)
	readr::write_tsv(x = signature_df, file = as.character(opt$output_file), col_names = TRUE, append = FALSE)
}
