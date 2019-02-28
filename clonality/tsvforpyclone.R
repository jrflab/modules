#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

optList = list(make_option("--sample_set", default = NULL, help = "sample set name"),
			   make_option("--normal_samples", default = NULL, help = "normal sample names"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options


tumor_samples = unlist(strsplit(opt$sample_set, split="_", fixed=TRUE))
normal_sample = unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE))
normal_sample = tumor_samples[tumor_samples %in% normal_sample]
tumor_samples = tumor_samples[!(tumor_samples %in% normal_sample)]

mutation_summary = read_tsv(file=paste0("sufam/", opt$sample_set, ".tsv"))
for (i in 1:length(tumor_samples)) {
	mutation_id = paste0(mutation_summary$Gene_Symbol, "_", mutation_summary$HGVSp)
	fsq = mutation_summary %>%
		  .[[paste0("MAF_", tumor_samples[i])]]
	qt = mutation_summary %>%
		  .[[paste0("qt_", tumor_samples[i])]]
	n = mutation_summary %>%
		.[[paste0("DP_", tumor_samples[i])]]
	flag = mutation_summary %>%
		   .[[paste0("CALL_", tumor_samples[i])]]
	
	fsq[flag==0] = 0
	var_counts = round(fsq*n)
	ref_counts = round((1-fsq)*n)
	normal_cn = rep(2, length(mutation_id))
	minor_cn = rep(0, length(mutation_id))
	major_cn = qt
	sample_summary = data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
	write.table(sample_summary, paste0("pyclone/", opt$sample_set, "/", tumor_samples[i], ".tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
}
