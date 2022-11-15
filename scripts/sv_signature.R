#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_names", default = NA, type = 'character', help = "sample names"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"),
	       make_option("--p_value", default = "0.05", type = 'character', help = "cluster sv p-value"),
	       make_option("--n_sv", default = "50", type = 'character', help = "number of sv"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_names = as.character(opt$sample_names)
	bedpe_org = readr::read_tsv(file = paste0("sv_signature/", sample_names, "/", sample_names, ".merged.bedpe"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert()
	bedpe_cli = readr::read_tsv(file = paste0("sv_signature/", sample_names, "/", sample_names, ".sv_clusters_and_footprints.tsv"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert() %>%
		    dplyr::select(chrom1 = X1,
				  start1 = X2,
				  end1 = X3,
				  chrom2 = X4,
				  start2 = X5,
				  end2 = X6,
				  sv_id = X7,
				  pe_support = X8,
				  strand1 = X9,
				  strand2 = X10,
				  n_svs = X12,
				  p_value = X17) %>%
		    dplyr::mutate(p_value = case_when(
			    is.na(p_value) ~ 1,
			    TRUE ~ p_value
		    )) %>%
		    dplyr::left_join(bedpe_org %>%
				     dplyr::select(chrom1, start1, end1, chrom2, start2, end2, svclass),
				     by = c("chrom1", "start1", "end1", "chrom2", "start2", "end2")) %>%
		    dplyr::filter(!is.na(svclass)) %>%
		    dplyr::mutate(is_clustered = case_when(
			    p_value < as.numeric(opt$p_value) & n_svs > 2*as.numeric(opt$n_sv) ~ "Cplx2",
			    p_value < as.numeric(opt$p_value) & n_svs > as.numeric(opt$n_sv) ~ "Cplx1",
			    TRUE ~ ""
		    )) %>%
		    dplyr::mutate(svclass = case_when(
			    svclass == "TRA" & is_clustered != "" ~ paste0(is_clustered, svclass),
			    TRUE ~ svclass
		    )) %>%
		    dplyr::select(chrom1, start1, end1, chrom2, start2, end2, sv_id, pe_support, strand1, strand2, svclass)
	
	write_tsv(x = bedpe_cli, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==2) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_names), split = " ", fixed=TRUE))
	feature_counts = list()
	for (i in 1:length(sample_names)) {
		feature_counts[[i]] = readr::read_tsv(file = paste0("sv_signature/", sample_names[i], "/", sample_names[i], ".merged.txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
				      readr::type_convert() %>%
				      dplyr::rename(sv_class = X1,
					            sv_count = manual_sv_type) %>%
				      dplyr::mutate(sample_name = sample_names[i])
	}
	feature_counts = do.call(bind_rows, feature_counts)
	write_tsv(x = feature_counts, path = as.character(opt$output_file), append = FALSE, col_names = TRUE)
	
}
