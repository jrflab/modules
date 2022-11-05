#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("superheat"))
suppressPackageStartupMessages(library("RColorBrewer"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list(make_option("--option", default = NA, type = 'character', help = "analysis type"),
               make_option("--sample_set", default = NA, type = 'character', help = "sample set"),
	       make_option("--normal_sample", default = NA, type = 'character', help = "normal sample"),
	       make_option("--input_file", default = NA, type = 'character', help = "input file"),
	       make_option("--output_file", default = NA, type = 'character', help = "output file"),
	       make_option("--num_iter", default = NA, type = 'character', help = "mcmc iterations"))
parser = OptionParser(usage = "%prog", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$option) == 1) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = "_"))
	normal_sample = as.character(opt$normal_sample)
	sample_set = setdiff(sample_set, normal_sample)
	pyclone = list()
	for (i in 1:length(sample_set)) {
		sufam = readr::read_tsv(file = paste0("pyclone_13/", sample_set[i], "/", sample_set[i], ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
				      var_counts = t_alt_count,
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
			       dplyr::select(mutation_id, sample_id, ref_counts, var_counts, normal_cn, major_cn, minor_cn)
		
	}
	pyclone = do.call(rbind, pyclone) %>%
		  dplyr::filter(!is.na(ref_counts)) %>%
		  dplyr::filter(!is.na(var_counts)) %>%
		  dplyr::mutate(var_counts = ifelse(var_counts<=3, 0, var_counts)) %>%
		  dplyr::filter((ref_counts+var_counts)>10) %>%
		  dplyr::filter(!is.na(major_cn)) %>%
		  dplyr::filter(major_cn != 0) %>%
		  dplyr::mutate(minor_cn = ifelse(is.na(minor_cn), 0, minor_cn))
	
	smry = pyclone %>%
	       dplyr::group_by(mutation_id) %>%
	       dplyr::summarize(n = n()) %>%
	       dplyr::ungroup()
	
	pyclone = pyclone %>%
		  dplyr::left_join(smry, by = "mutation_id") %>%
		  dplyr::filter(n == length(sample_set))
	
	for (i in 1:length(sample_set)) {
		pyclone_ft = pyclone %>%
			     dplyr::filter(sample_id == sample_set[i]) %>%
			     dplyr::select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
		readr::write_tsv(x = pyclone_ft, file = paste0("pyclone_13/", opt$sample_set, "/", sample_set[i], ".tsv"), append = FALSE, col_names = TRUE)
	}
	
} else if (as.numeric(opt$option) == 2) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = "_"))
	normal_sample = as.character(opt$normal_sample)
	sample_set = setdiff(sample_set, normal_sample)
	params = list()
	for (i in 1:length(sample_set)) {
		params[[i]] = readr::read_tsv(file = paste0("facets/cncf/", sample_set[i], "_", normal_sample, ".out"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			      readr::type_convert() %>%
			      dplyr::filter(grepl("# Purity", X1)) %>%
			      dplyr::mutate(X1 = gsub("# Purity = ", "", X1)) %>%
			      readr::type_convert() %>%
			      .[["X1"]]
	}
	cat(paste0("num_iters: ", as.numeric(opt$num_iter), "\n\n"), file = as.character(opt$output_file), append = FALSE)
	cat("base_measure_params:\n", file = as.character(opt$output_file), append = TRUE)
	cat("  alpha: 1\n", file = as.character(opt$output_file), append = TRUE)
	cat("  beta: 1\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("concentration:\n", file = as.character(opt$output_file), append = TRUE)
	cat("  value: 1.0\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("  prior:\n", file = as.character(opt$output_file), append = TRUE)
	cat("    shape: 1.0\n", file = as.character(opt$output_file), append = TRUE)
	cat("    rate: 0.001\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
    	cat("density: pyclone_beta_binomial\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("beta_binomial_precision_params:\n", file = as.character(opt$output_file), append = TRUE)
	cat("  value: 1000\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("  prior:\n", file = as.character(opt$output_file), append = TRUE)
  	cat("    shape: 1.0\n", file = as.character(opt$output_file), append = TRUE)
	cat("    rate: 0.0001\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("  proposal:\n", file = as.character(opt$output_file), append = TRUE)
	cat("    precision: 0.01\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("working_dir: pyclone_13/", file = as.character(opt$output_file), append = TRUE)
	cat(as.character(opt$sample_set), file = as.character(opt$output_file), append = TRUE)
	cat("\n\n", file = as.character(opt$output_file), append = TRUE)
	cat("trace_dir: trace\n", file = as.character(opt$output_file), append = TRUE)
	cat("init_method: connected\n", file = as.character(opt$output_file), append = TRUE)
	cat("\n", file = as.character(opt$output_file), append = TRUE)
	cat("samples:\n", file = as.character(opt$output_file), append = TRUE)
	for (i in 1:length(sample_set)) {
		cat(paste0("  ", sample_set[i], ":\n"), file = as.character(opt$output_file), append = TRUE)
		cat(paste0("    mutations_file: pyclone_13/", as.character(opt$sample_set), "/", sample_set[i], ".yaml\n\n"), file = as.character(opt$output_file), append = TRUE)
		cat("    tumour_content:\n", file = as.character(opt$output_file), append = TRUE)
		cat(paste0("      value: ", params[[i]], "\n"), file = as.character(opt$output_file), append = TRUE)
		cat("\n", file = as.character(opt$output_file), append = TRUE)
		cat("    error_rate: 0.01\n", file = as.character(opt$output_file), append = TRUE)
		cat("\n", file = as.character(opt$output_file), append = TRUE)
	}
} else if (as.numeric(opt$option) == 3) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " "))
	pyclone = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		  readr::type_convert()
	
	pyclone_ft = list()
	index = 1
	for (i in 1:(length(sample_set)-1)) {
		for (j in (i+1):length(sample_set)) {
			pyclone_ft[[index]] = pyclone %>%
					      dplyr::filter(sample_id == sample_set[i]) %>%
					      dplyr::select(mutation_id,
							    cluster_id,
							    sample_id_x = sample_id,
							    cellular_prevalence_x = cellular_prevalence,
							    cellular_prevalence_std_x = cellular_prevalence_std) %>%
					      dplyr::full_join(pyclone %>%
							       dplyr::filter(sample_id == sample_set[j]) %>%
							       dplyr::select(mutation_id,
									     sample_id_y = sample_id,
									     cellular_prevalence_y = cellular_prevalence,
									     cellular_prevalence_std_y = cellular_prevalence_std),
							       by = "mutation_id") %>%
					      readr::type_convert()
			index = index + 1
		}
	}
	pyclone_ft = do.call(bind_rows, pyclone_ft) %>%
		     readr::type_convert()
	
	smry_c = pyclone_ft %>%
		 dplyr::group_by(mutation_id) %>%
		 dplyr::summarize(cluster_id = unique(cluster_id)) %>%
		 dplyr::ungroup() %>%
		 dplyr::group_by(cluster_id) %>%
		 dplyr::summarize(n = n())
	smry_p = pyclone %>%
		 dplyr::group_by(cluster_id, sample_id) %>%
		 dplyr::summarize(mean_cellular_prevalence = mean(cellular_prevalence)) %>%
		 dplyr::ungroup() %>%
		 dplyr::group_by(cluster_id) %>%
		 dplyr::summarize(is_clonal = max(mean_cellular_prevalence))
	
	pyclone_ft = pyclone_ft %>%
		     dplyr::left_join(smry_c, by = "cluster_id") %>%
		     dplyr::left_join(smry_p, by = "cluster_id")
	
	smry_cl = pyclone %>%
		  dplyr::group_by(cluster_id) %>%
		  dplyr::summarize(mean = mean(cellular_prevalence)) %>%
		  dplyr::ungroup() %>%
		  dplyr::arrange(desc(mean)) %>%
		  dplyr::mutate(cluster_id_ordered = nrow(.):1)
	
	pyclone_ft = pyclone_ft %>%
		     dplyr::left_join(smry_cl, by = "cluster_id")
	
	colourCount = length(unique(pyclone_ft$cluster_id_ordered))
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
		
	plot_ = pyclone_ft %>%
		ggplot(aes(x = 100*cellular_prevalence_x, y = 100*cellular_prevalence_y, color = factor(cluster_id_ordered), size = n)) +
		geom_point(stat = "identity", alpha = .75, shape = 21) +
		scale_color_manual(values = getPalette(colourCount)) +
		xlab("\n\nCCF (%)\n") +
		ylab("\nCCF (%)\n\n") +
		guides(color = guide_legend(title = "Cluster"),
		       size = guide_legend(title = "N")) +
		facet_wrap(~sample_id_x+sample_id_y)
	pdf(file = as.character(opt$output_file), width = 21, height = 21)
	print(plot_)
	dev.off()
	
} else if (as.numeric(opt$option) == 4) {
	sample_set = unlist(strsplit(x = as.character(opt$sample_set), split = " "))
	pyclone = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		  readr::type_convert() %>%
		  dplyr::mutate(sample_id = paste0(sample_id, "     "))
	
	pyclone_mt = pyclone %>%
		     reshape2::dcast(formula = mutation_id~sample_id, value.var = "cellular_prevalence") %>%
		     dplyr::left_join(pyclone %>%
				      dplyr::select(mutation_id, cluster_id) %>%
				      dplyr::filter(!duplicated(mutation_id)), by = "mutation_id")
	
	smry_cl = pyclone %>%
		  dplyr::group_by(cluster_id) %>%
		  dplyr::summarize(mean = mean(cellular_prevalence)) %>%
		  dplyr::ungroup() %>%
		  dplyr::arrange(desc(mean)) %>%
		  dplyr::mutate(cluster_id_ordered = nrow(.):1)
	
	pyclone_mt = pyclone_mt %>%
		     dplyr::left_join(smry_cl, by = "cluster_id")
	
	index = order(apply(pyclone_mt %>% dplyr::select(-mutation_id, -cluster_id, -cluster_id_ordered), 1, mean), decreasing = TRUE)
	pyclone_mt = pyclone_mt[index,,drop=FALSE]
	pyclone_mt = pyclone_mt %>%
		     dplyr::arrange(cluster_id_ordered)
		
	
	pdf(file = as.character(opt$output_file), width = 10, height = 21)
	superheat(X = pyclone_mt %>%
		      dplyr::select(-mutation_id, -cluster_id, -cluster_id_ordered, -mean),
		  membership.rows = pyclone_mt %>% .[["cluster_id_ordered"]],
		  pretty.order.rows = FALSE,
		  pretty.order.cols = TRUE,
		  row.dendrogram = FALSE,
		  col.dendrogram = FALSE,
		  smooth.heat = FALSE,
		  scale = FALSE,
		  heat.pal = c(rep("#d9d9d9", 10), "#9ecae1", "#4292c6", "#2171b5", "#08519c", "#08306b"),
		  legend = FALSE,
		  grid.hline = FALSE,
		  grid.vline = TRUE,
		  force.grid.hline = TRUE,
		  force.grid.vline = TRUE,
		  grid.hline.col = "white",
		  grid.vline.col = "white",
		  grid.hline.size = .05,
		  grid.vline.size = 1,
		  bottom.label.text.angle = 90,
		  bottom.label.text.alignment = "right")
	dev.off()

}
