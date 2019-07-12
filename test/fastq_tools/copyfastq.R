#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list (
				make_option("--sample_name", default = NULL, help = "sample name"),
				make_option("--fastq_files", default = NULL, help = "fastq files")
			   )

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

sample_name = opt$sample_name
fastq_files = unlist(strsplit(opt$fastq_files, split=" ", fixed=TRUE))

file.copy(from=fastq_files[1], to=paste0("marianas/", sample_name, "/", sample_name, "_R1.fastq.gz"))
file.copy(from=fastq_files[2], to=paste0("marianas/", sample_name, "/", sample_name, "_R2.fastq.gz"))
