#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file_in", default = NA, type = 'character', help = "input file name"),
				  make_option("--file_out", default = NA, type = 'character', help = "output file name"))

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load(opt$file_in)
alpha = ifelse(is.na(fit$purity), 1, fit$purity)
psi = ifelse(is.na(fit$ploidy), 2, fit$ploidy)
gamma = 1
x = fit$cncf[,"cnlr.median"]
absolute_copies = round(((((2^(x/gamma))*(alpha*psi+(1-alpha)*2)) - ((1-alpha)*2))/alpha))
index = absolute_copies!=round(psi)
if (sum(index, na.rm=TRUE)!=0) {
	genome_footprint = sum(as.numeric(fit$cncf[,"end"]-fit$cncf[,"start"]), na.rm=TRUE)
	genome_altered = sum(as.numeric(fit$cncf[index,"end"]-fit$cncf[index,"start"]), na.rm=TRUE)/genome_footprint
} else {
	genome_altered = 0
}
cat(paste0(gsub("facets/cncf/","", gsub(".Rdata", "", opt$file_in)), "\t", genome_altered), file = opt$file_out, append=FALSE)
cat("\n", file = opt$file_out, append=TRUE)

warnings()
