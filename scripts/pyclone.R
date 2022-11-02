#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))

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

if (as.numeric(opt$option) == 1) {
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
	
} else if (as.numeric(opt$option) == 2) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " "))
	pyclone = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		  readr::type_convert()
	
	pyclone_ft = list()
	index = 1
	for (i in 1:(length(sample_set)-1)) {
		for (j in 2:length(sample_set)) {
			pyclone_ft[[index]] = pyclone %>%
					      dplyr::filter(sample_id == sample_set[i]) %>%
					      dplyr::rename(sample_id_x = sample_id,
							    cellular_prevalence_x = cellular_prevalence,
							    cellular_prevalence_std_x = cellular_prevalence_std) %>%
					      dplyr::full_join(pyclone %>%
							       dplyr::filter(sample_id == sample_set[j]) %>%
							       dplyr::rename(sample_id_y = sample_id,
									     cellular_prevalence_y = cellular_prevalence,
									     cellular_prevalence_std_y = cellular_prevalence_std)) %>%
					      readr::type_convert()
			index = index + 1
		}
	}
	pyclone_ft = do.call(bind_rows, pyclone_ft) %>%
		     dplyr::filter(cellular_prevalence_x > 0 & cellular_prevalence_y > 0)
	smry_c = pyclone_ft %>%
		 dplyr::group_by(mutation_id) %>%
		 dplyr::summarize(cluster_id = unique(cluster_id)) %>%
		 dplyr::ungroup() %>%
		 dplyr::group_by(cluster_id) %>%
		 dplyr::summarize(n = n())
	pyclone_ft = pyclone_ft %>%
		     dplyr::left_join(smry_c, by = "cluster_id")
	
	plot_ = pyclone_ft %>%
		ggplot(aes(x = 100*cellular_prevalence_x, y = 100*cellular_prevalence_y, color = factor(cluster_id), size = n)) +
		geom_point(stat = "identity", alpha = .75, shape = 21) +
		scale_color_brewer(type = "qual", palette = 6) +
		xlab("\n\nCCF (%)\n") +
		ylab("\nCCF (%)\n\n") +
		guides(color = guide_legend(title = "Cluster"),
		       size = guide_legend(title = "N")) +
		facet_wrap(sample_id_x~sample_id_y)
	pdf(file = as.character(opt$output_file), width = 21, height = 21)
	print(plot_)
	dev.off()
		
		
	
	
}
