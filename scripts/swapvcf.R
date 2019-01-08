#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file", default = NA, type = 'character', help = "input file name"),
				  make_option("--tumor", default = NA, type = 'character', help = "tumor sample name"),
				  make_option("--normal", default = NA, type = 'character', help = "normal sample name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

vcf = read.table(file=opt$file, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
n = ncol(vcf)
n1 = vcf[,n,drop=TRUE]
n2 = vcf[,n-1,drop=TRUE]
vcf[,n-1] = n1
vcf[,n] = n2
colnames(vcf) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", opt$tumor, opt$normal)
system(paste0("grep '##' ", opt$file, " > ", opt$file, ".tmp"))
write.table(vcf, file=paste0(opt$file, ".tmp"), append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
warnings()
