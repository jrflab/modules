#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file name"),
				  make_option("--option", default = NA, type = 'character', help = "0/1 on- or offtarget"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

in_file = unlist(strsplit(x=opt$in_file, split=" ", fixed=TRUE))
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

'IQR' <- function(x)
{
	x = na.omit(x)
	iq = stats::IQR(abs(x[1:(length(x)-1)] - x[2:length(x)]))
	return(invisible(iq))
}

data = matrix(NA, nrow=length(in_file), ncol=3)
for (i in 1:length(in_file)) {
	print(i)
	data = read.csv(file=in_file[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% 1:22 & data[,"gene"] == ifelse(opt$option==1, "-", "Antitarget")
	data[i,1] = MAD(data[index,"log2"])
	data[i,2] = MAPD(data[index,"log2"])
	data[i,3] = IQR(data[index,"log2"])
}
colnames(data) = c("MAD", "MAPD", "IQR")
write.table(data, file=opt_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
