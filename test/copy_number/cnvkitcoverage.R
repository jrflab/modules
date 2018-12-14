#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))

'fix' <- function(x)
{
	y = log2(x)
	m = Mclust(x, G=1:2)
	for (i in 1:m$G) {
		index = m$classification==i
		z = y[index]
		z = (z - mean(z, na.rm=TRUE))/sd(z, na.rm=TRUE)
		y[index] = z
	}
	y = 2^y
	return(invisible(y))
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--file", default = NA, type = 'character', help = "file name and path"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

filein = opt$file
fileout = gsub(pattern=".tmp", replacement=".cnn", x=filein, fixed=TRUE)

data = read.csv(file=filein, header=TRUE, sep="\t", stringsAsFactors=FALSE)
index = data[,1]=="X" | data[,1]=="Y"
x = fix(data[!index,"depth"])
y = log2(x)
inx = x<0 | is.infinite(y) | is.na(x) | is.na(y)
x[inx] = 0
y[inx] = -20
data[!index,"depth"] = x
data[!index,"log2"] = y
write.table(data, file=fileout, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
