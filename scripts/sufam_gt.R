#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

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


all_vars = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
CHROM = all_vars[,"CHROM"]
POS = all_vars[,"POS"]
ID = all_vars[,"ID"]
REF = all_vars[,"REF"]
ALT = all_vars[,"ALT"]
QUAL = FILTER = rep(".", nrow(all_vars))
INFO = paste0(all_vars[,"SYMBOL"], all_vars[,"HGVSp_Short"])
vcf = data.frame(CHROM, POS, ID, REF, ALT, QUAL, INFO)

cat("#", file="sufam/pdx.vcf", append=FALSE)
write.table(vcf, file="sufam/pdx.vcf", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
