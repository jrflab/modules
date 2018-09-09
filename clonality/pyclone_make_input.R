#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--file_name", default = NULL, help = "sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options
mutation_summary = read.csv(file=opt$file_name, header=TRUE, sep="\t", stringsAsFactors=FALSE)
index = grep("DP", colnames(mutation_summary))[-1]
sample_names = gsub("DP_", "", x=colnames(mutation_summary)[index], fixed=TRUE)

mutation_id = paste0(mutation_summary[,"Gene_Symbol"], "_", mutation_summary[,"HGVSp"])

for (i in 1:length(sample_names)) {
	ref_counts = round((1-mutation_summary[,paste0("MAF_", sample_names[i])])*mutation_summary[,paste0("DP_", sample_names[i])])
	var_counts = round((mutation_summary[,paste0("MAF_", sample_names[i])])*mutation_summary[,paste0("DP_", sample_names[i])])
	normal_cn = rep(2, length(mutation_id))
	minor_cn = mutation_summary[,paste0("qt_", sample_names[i])] - mutation_summary[,paste0("q2_", sample_names[i])]
	major_cn = mutation_summary[,paste0("q2_", sample_names[i])]
	sample_summary = data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
	write.table(sample_summary, paste0(gsub(".tsv", "/", gsub("sufam/", "pyclone/", opt$file_name)), sample_names, ".tsv"))
}

