#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--normal_files", default = NA, type = 'character', help = "normal input files"),
				  make_option("--tumor_files", default = NA, type = 'character', help = "tumor input files"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file"),
				  make_option("--option", default = NA, type = 'character', help = "1-0 for ontarget or offtarget"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

normal_files = unlist(strsplit(x=opt$normal_files, split=" ", fixed=TRUE))
normal_samples = gsub(pattern=".cnr", replacement="", x=gsub(pattern="cnvkit/cnr/", replacement="", x=normal_files, fixed=TRUE), fixed=TRUE)
tumor_files = unlist(strsplit(x=opt$tumor_files, split=" ", fixed=TRUE))
tumor_samples = gsub(pattern=".cnr", replacement="", x=gsub(pattern="cnvkit/cnr/", replacement="", x=tumor_files, fixed=TRUE), fixed=TRUE)
out_file = opt$out_file

'MAD' <- function(x)
{
	x = na.omit(x)
	q2 = mad(x)
	return(invisible(q2))
}

'MAPD' <- function(x)
{
	x = na.omit(x)
	q2 = median(abs(x[1:(length(x)-1)] - x[2:length(x)]))
	return(invisible(q2))
}

'MIQR' <- function(x)
{
	x = na.omit(x)
	iq = stats::IQR(abs(x[1:(length(x)-1)] - x[2:length(x)]))
	return(invisible(iq))
}

data = matrix(NA, nrow=length(c(normal_samples, tumor_samples)), ncol=3, dimnames=list(c(normal_samples, tumor_samples), c("MAD", "MAPD", "IQR")))
for (i in 1:length(normal_files)) {
	print(i)
	data = read.csv(file=normal_files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% 1:22 & data[,"gene"] == ifelse(opt$option==1, "-", "Antitarget")
	data[normal_samples[i],1] = MAD(data[index,"log2"])
	data[normal_samples[i],2] = MAPD(data[index,"log2"])
	data[normal_samples[i],3] = MIQR(data[index,"log2"])
}
for (i in 1:length(tumor_files)) {
	print(i)
	data = read.csv(file=tumor_files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% 1:22 & data[,"gene"] == ifelse(opt$option==1, "-", "Antitarget")
	data[tumor_samples[i],1] = MAD(data[index,"log2"])
	data[tumor_samples[i],2] = MAPD(data[index,"log2"])
	data[tumor_samples[i],3] = MIQR(data[index,"log2"])
}

colnames(data) = c("MAD", "MAPD", "IQR")
#data = cbind("SAMPLE_NAME"=c(normal_samples, tumor_samples), "SAMPLE_TYPE"=c(rep("N", length(normal_samples)), rep("T", length(tumor_samples))), data)
write.table(data, file=out_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
