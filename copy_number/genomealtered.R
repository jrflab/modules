#!/usr/bin/env Rscript
# compute the fraction of genome altered

suppressPackageStartupMessages(library("optparse"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file_in", default = NA, type = 'character', help = "input file name"),
				  make_option("--file_out", default = NA, type = 'character', help = "output file name"))

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load(opt$file_in)
alpha = fit$purity
psi = fit$ploidy
gamma = 1
x = fit$cncf[,"cnlr.median"]
absolute_copies = round(((((2^(x/gamma))*(alpha*psi+(1-alpha)*2)) - ((1-alpha)*2))/alpha))
index = absolute_copies!=round(psi)
genome_footprint = sum(fit$cncf[!is.na(index),"end"]-fit$cncf[!is.na(index),"start"], na.rm=TRUE)
genome_altered = sum((fit$cncf[index,"end"] - fit$cncf[index,"start"]), na.rm=TRUE)/genome_footprint
genome_altered = ifelse(is.nan(genome_altered), 0, genome_altered)
cat(paste0(gsub("facets/cncf/","", gsub(".Rdata", "", opt$file_in)), "\t", genome_altered), file = opt$file_out, append=FALSE)
cat("\n", file = opt$file_out, append=TRUE)

warnings()