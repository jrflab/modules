#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("Starfish"))


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
	bed = readr::read_tsv(file = as.character(opt$input_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::rename(chrom1 = X1,
			    start1 = X2,
			    end1 = X3,
			    chrom2 = X4,
			    start2 = X5,
			    end2 = X6,
			    sv_id = X7,
			    pe_support = X8,
			    strand1 = X9,
			    strand2 = X10,
			    svclass = X11) %>%
	      dplyr::select(chrom1, pos1 = start1, chrom2, pos2 = start2, strand1, strand2, svtype = svclass) %>%
	      dplyr::mutate(svtype = case_when(
		      svtype == "INV" & strand1 == "+" & strand2 == "+" ~ "h2hINV",
		      svtype == "INV" & strand1 == "-" & strand2 == "-" ~ "t2tINV",
		      TRUE ~ svtype
	      )) %>%
	      dplyr::mutate(sample = sample_name)
	readr::write_tsv(x = bed, file = as.character(opt$output_file), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==2) {
	sample_name = as.character(opt$sample_name)
	data = readr::read_tsv(file = as.character(opt$input_file), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       dplyr::select(chromosome = chrom,
			     start = loc.start,
			     end = loc.end,
			     total_cn = tcn.em) %>%
	       dplyr::mutate(sample = sample_name)
	readr::write_tsv(x = data, file = as.character(opt$output_file), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option)==3) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_name), split = " ", fixed = TRUE))
	sv_df = cn_df = gd_df = list()
	for (i in 1:length(sample_names)) {
		sv_df[[i]] = readr::read_tsv(file = paste0("star_fish/", sample_names[i], "/", sample_names[i], ".merged_sv.bedpe"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
			     readr::type_convert()
		cn_df[[i]] = readr::read_tsv(file = paste0("star_fish/", sample_names[i], "/", sample_names[i], ".merged_cn.txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
		     	     readr::type_convert()
		gd_df[[i]] = dplyr::tibble(sample = sample_names[i], gender = "unknown") %>%
		     	     readr::type_convert()
	}
	sv_df = do.call(rbind, sv_df)
	cn_df = do.call(rbind, cn_df)
	gd_df = do.call(rbind, gd_df)
	starfish_link_out = starfish_link(sv_file = sv_df, prefix = "star_fish/summary/")
	if (length(starfish_link_out) == 1) {
		cat(starfish_link_out, file = as.character(opt$output_file), append = FALSE)
	} else {
		starfish_feature_out = starfish_feature(cgr = starfish_link_out$starfish_call,
							complex_sv = starfish_link_out$interleave_tra_complex_sv,
							cnv_file = cn_df %>%
								   dplyr::mutate(chromosome = as.character(chromosome)) %>%
								   dplyr::mutate(chromosome = case_when(
									   chromosome == "23" ~ "X",
									   TRUE ~ chromosome)),
							gender_file = gd_df,
							prefix = "star_fish/summary/",
							genome_v = "hg19",
							cnv_factor = "auto",
							arm_del_rm = TRUE)
		starfish_sig_out = starfish_sig(cluster_feature = starfish_feature_out$cluster_feature,
					        prefix = "star_fish/summary/",
					        cmethod = "class")
		wd = getwd()
		setwd("star_fish/summary/")
		starfish_plot(sv_file = sv_df, cnv_file = cn_df, cgr = starfish_link_out$starfish_call, genome_v = "hg19")
		setwd(wd)
		cat("taskcomplete!!", file = as.character(opt$output_file), append = FALSE)
	}

} else if (as.numeric(opt$option)==4) {
	df = readr::read_csv(file = "star_fish/summary/_pcawg_6signatures_class.csv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
	     dplyr::mutate(sample_name = unlist(lapply(cluster_id, function(x) { paste0(unlist(strsplit(x, "_", fixed = TRUE))[1:2], collapse="_")})))
	readr::write_tsv(x = df, file = as.character(opt$output_file), append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$option)==5) {
	df = readr::read_csv(file = "star_fish/summary/_CGR_feature_matrix.csv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
	     dplyr::rename(sample_name = sample)
	readr::write_tsv(x = df, file = as.character(opt$output_file), append = FALSE, col_names = TRUE)
}
