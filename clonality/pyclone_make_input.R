#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--file_name", default = NULL, help = "file name"),
			   make_option("--sample_name", default = NULL, help = "sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options
mutation_summary = read.csv(file=opt$file_name, header=TRUE, sep="\t", stringsAsFactors=FALSE)
mutation_id = paste0(mutation_summary[,"Gene_Symbol"], "_", mutation_summary[,"HGVSp"])

ref_counts = round((1-mutation_summary[,paste0("MAF_", opt$sample_name)])*mutation_summary[,paste0("DP_", opt$sample_name)])
var_counts = round((mutation_summary[,paste0("MAF_", opt$sample_name)])*mutation_summary[,paste0("DP_", opt$sample_name)])
normal_cn = rep(2, length(mutation_id))
minor_cn = mutation_summary[,paste0("qt_", opt$sample_name)] - mutation_summary[,paste0("q2_", opt$sample_name)]
major_cn = mutation_summary[,paste0("q2_", opt$sample_name)]
sample_summary = data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
write.table(sample_summary, paste0(gsub(".tsv", "/", gsub("sufam/", "pyclone/", opt$file_name)), opt$sample_name, ".tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)

