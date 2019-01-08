#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

in_file_names = unlist(strsplit(x=opt$in_file, split=" ", fixed=TRUE))
sample_names = gsub(".txt", "", x=gsub("sufam/", "", in_file_names, fixed=TRUE), fixed=TRUE)
out_file_name = opt$out_file
DP = AD = MAF = list()
for (i in 1:length(in_file_names)) {
	tmp = read.csv(file=in_file_names[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	DP[[i]] = tmp[,"cov"]
	AD[[i]] = tmp[,"val_al_count"]
	MAF[[i]] = tmp[,"val_maf"]
}
DP = do.call(cbind, DP)
colnames(DP) = paste0("DP_", sample_names)
AD = do.call(cbind, AD)
colnames(AD) = paste0("AD_", sample_names)
MAF = do.call(cbind, MAF)
colnames(MAF) = paste0("MAF_", sample_names)
vcf = read.table(file="sufam/PDX.vcf", header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
chr = vcf[,1]
pos = vcf[,2]
ref = vcf[,4]
alt = vcf[,5]
gene_symbol = unlist(lapply(strsplit(vcf[,7], "p.", fixed=TRUE), function(x) { x[1] }))
hgvsp_short = paste0("p.", unlist(lapply(strsplit(vcf[,7], "p.", fixed=TRUE), function(x) { x[2] })))
res = cbind(chr, pos, ref, alt, gene_symbol, hgvsp_short, DP, AD, MAF)
colnames(res)[1:6] = c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele", "Gene_Symbol", "HGVSp")
write.table(res, file=out_file_name, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
