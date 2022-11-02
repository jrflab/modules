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
	       make_option("--normal_sample", default = NA, type = 'character', help = "normal sample"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option)==1) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = "_"))
	normal_sample = as.character(opt$normal_sample)
	sample_set = setdiff(sample_set, normal_sample)
	pyclone = list()
	for (i in 1:length(sample_set)) {
		sufam = readr::read_tsv(file = paste0("pyclone/", sample_set[i], ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(Chromosome = chrom,
				      Position = pos,
				      Reference_Allele = val_ref,
				      Alternate_Allele = val_alt,
				      t_depth = cov,
				      t_alt_count = val_al_count) %>%
		 	dplyr::mutate(t_ref_count = t_depth - t_alt_count) %>%
			dplyr::mutate(mutation_id = paste0(Chromosome, ":", Position, ":", Reference_Allele, ":", Alternate_Allele),
				      ref_counts = t_ref_count,
				      alt_counts = t_alt_count,
				      normal_cn = 2)
		
		facets = readr::read_tsv(file = paste0("facets/cncf/", sample_set[i], "_", normal_sample, ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 dplyr::mutate(chrom = ifelse(chrom == "23", "X", chrom)) %>%
			 dplyr::mutate(Chromosome = chrom,
				       Start_Position = loc.start,
				       End_Position = loc.end,
				       minor_cn = ifelse(is.na(lcn.em), "0", lcn.em),
				       major_cn = tcn.em) %>%
		 readr::type_convert() %>%
		 dplyr::mutate(major_cn = major_cn - minor_cn) %>%
		 dplyr::select(Chromosome, Start_Position, End_Position, minor_cn, major_cn)
		 
		pyclone[[i]] = sufam %>%
		  	       dplyr::mutate(Chromosome = ifelse(Chromosome == "X", "23", Chromosome)) %>%
			       dplyr::mutate(Start_Position = Position,
					     End_Position = Position +1) %>%
			       readr::type_convert() %>%
			       fuzzyjoin::genome_left_join(facets %>%
							   dplyr::mutate(Chromosome = ifelse(Chromosome == "X", "23", Chromosome)) %>%
							   readr::type_convert(),
							   by = c("Chromosome", "Start_Position", "End_Position")) %>%
			       dplyr::mutate(sample_id = sample_set[i]) %>%
			       dplyr::select(mutation_id, sample_id, ref_counts, alt_counts, normal_cn, major_cn, minor_cn)
	
		params = readr::read_tsv(file = paste0("facets/cncf/", sample_set[i], "_", normal_sample, ".out"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::filter(grepl("# Purity", X1)) %>%
			 dplyr::mutate(X1 = gsub("# Purity = ", "", X1)) %>%
			 readr::type_convert() %>%
			 .[["X1"]]
		
		pyclone[[i]] = pyclone[[i]] %>%
			       dplyr::mutate(tumour_content = params)
	}
	pyclone = do.call(rbind, pyclone) %>%
		  dplyr::filter(!is.na(ref_counts)) %>%
		  dplyr::filter(!is.na(alt_counts)) %>%
		  dplyr::filter(!is.na(major_cn)) %>%
		  dplyr::filter(!is.na(minor_cn))
	
	readr::write_tsv(x = pyclone, file = opt$output_file, append = FALSE, col_names = TRUE)
}
