#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))

optList = list(make_option("--sample_name", default = NULL, help = "tumor normal sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_sample = unlist(strsplit(opt$sample_name, split="_", fixed=TRUE))[1]
normal_sample = unlist(strsplit(opt$sample_name, split="_", fixed=TRUE))[2]

file_paths = list(
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/pyclone.tsv"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/report.tsv")
)

mutation_summary = read_tsv(file=file_paths[[1]], col_types = cols(.default = col_character()), col_names = c("ID", "CCF", "STD", "C_ID")) %>%
				   type_convert() %>%
				   mutate(C_ID = factor(C_ID)) %>%
				   mutate(CCF = as.numeric(CCF)) %>%
				   mutate(STD = as.numeric(STD)) %>%
				   slice(-1)

tmp = mutation_summary %>%
	  group_by(C_ID) %>%
	  summarize(
	  		N = n(),
	  		Mean_CCF = mean(CCF),
	    	Median_CCF = median(CCF),
	    	Std_CCF = sd(CCF),
	    	Min_CCF = min(CCF),
	    	Max_CCF = max(CCF),
	    	Mean_Sd = mean(STD),
	    	Median_Sd = median(STD),
	    	Std_Sd = sd(STD),
	    	Min_Sd = min(STD),
	    	Max_Sd = max(STD)) %>%
	rename(Cluster_ID = C_ID)

write_tsv(x=tmp, path=file_paths[[2]])
