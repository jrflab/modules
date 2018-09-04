#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file_in", default = NA, type = 'character', help = "input file name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

vcf = read.csv(file=opt$file_in, header=FALSE, stringsAsFactors=FALSE)
n = ncol(vcf)
n1 = vcf[,n,drop=TRUE]
n2 = vcf[,n-1,drop=TRUE]
vcf[,n-1] = n1
vcf[,n] = n2
system(paste0("grep '#' ", opt$file_in, " > ", opt$file_in, ".tmp"))
write.table(vcf, file=paste0(opt$file_in, ".tmp"), append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
file.copy(from=paste0(opt$file_in, ".tmp"), to=opt$file_in, overwrite=TRUE)

warnings()
