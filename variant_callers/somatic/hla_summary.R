suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(make_option("--option", default = "NA", help = "which option?"),
                make_option("--sample_names", default = "NA", help = "sample names"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = unlist(strsplit(opt$sample_names, split=" ", fixed=TRUE))

if (as.numeric(opt$option)==1) {
	hla_genotypes = list()
	for (i in 1:length(sample_names)) {
		hla_genotypes[[i]] = readr::read_tsv(file = paste0("hla_polysolver/", sample_names[i], "/winners.hla.txt"),
						     col_names = FALSE, col_types = cols(.default = col_character())) %>%
				     readr::type_convert() %>%
				     dplyr::rename(hla = X1, major_allele = X2, minor_allele = X3) %>%
				     dplyr::mutate(sample_name = sample_names[i])
	}
	hla_genotypes = do.call(rbind, hla_genotypes)
	readr::write_tsv(x = hla_genotypes, path = "hla_polysolver/summary/hla_summary.txt", col_names = TRUE, append = FALSE)

} else if (as.numeric(opt$option)==2) {
	somatic_vars = list()
	for (i in 1:length(sample_names)) {
		somatic_vars[[i]] = readr::read_tsv(file = paste0("hla_polysolver/", sample_names[i], "/", sample_names[i], ".mutect.unfiltered.annotated"),
						    col_names = TRUE, col_types = cols(.default = col_character())) %>%
				    readr::type_convert()
	}
	somatic_vars = do.call(rbind, somatic_vars)
	if (nrow(somatic_vars)>0) {
		somatic_vars = somatic_vars %>%
			       dplyr::mutate(tumor_name = unlist(lapply(individual, function(x) { unlist(strsplit(x, split = "_", fixed = TRUE))[1] }))) %>%
			       dplyr::mutate(normal_name = unlist(lapply(individual, function(x) { unlist(strsplit(x, split = "_", fixed = TRUE))[2] })))
	}
	readr::write_tsv(x = somatic_vars, path = "hla_polysolver/summary/mutect_summary.txt", col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==3) {
	somatic_vars = list()
	for (i in 1:length(sample_names)) {
		somatic_vars[[i]] = readr::read_tsv(file = paste0("hla_polysolver/", sample_names[i], "/", sample_names[i], ".strelka_indels.unfiltered.annotated"),
						    col_names = TRUE, col_types = cols(.default = col_character())) %>%
				    readr::type_convert()
	}
	somatic_vars = do.call(rbind, somatic_vars)
	if (nrow(somatic_vars)>0) {
		somatic_vars = somatic_vars %>%
			       dplyr::mutate(tumor_name = unlist(lapply(individual, function(x) { unlist(strsplit(x, split = "_", fixed = TRUE))[1] }))) %>%
			       dplyr::mutate(normal_name = unlist(lapply(individual, function(x) { unlist(strsplit(x, split = "_", fixed = TRUE))[2] })))
	}
	readr::write_tsv(x = somatic_vars, path = "hla_polysolver/summary/strelka_summary.txt", col_names = TRUE, append = FALSE)

}
