#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("sleuth"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option('--annotation', type = 'character', default = NA, help = 'path to annotation file'),
               make_option('--samples', type = 'character', default = NA, help = 'list of samples names'))
parser = OptionParser(usage = "%prog",  option_list=optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

sample_names = unlist(strsplit(x=opt$samples, split=" ", fixed=TRUE))
annotation = readr::read_tsv(file=opt$annotation, col_names=TRUE, col_types=cols(.default=col_character()))
manifest = dplyr::tibble(sample = sample_names,
			 condition = rep(1, length(sample_names)),
			 path = paste0("kallisto/", sample_names))
data = sleuth::sleuth_prep(sample_to_covariates = manifest,
			   extra_bootstrap_summary = TRUE,
			   read_bootstrap_tpm = TRUE,
			   target_mapping = annotation,
			   aggregation_column = "hugo",
			   gene_mode = TRUE)
res = as.data.frame(sleuth_to_matrix(data, "obs_norm", "tpm"))
tpm_bygene = dplyr::tibble(gene_symbol = rownames(res)) %>%
	     dplyr::bind_cols(dplyr::as_tibble(res))
write_tsv(x=tpm_bygene, path="kallisto/tpm_by_gene.txt", append=FALSE, col_names=TRUE, quote_escape=FALSE)
