#!/usr/bin/env Rscript

for (lib in c("optparse", "dplyr")) {
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
}

optList <- list(make_option("--outFile", default = NULL, help = "output file"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (length(arguments$args) < 1) {
	cat("Need facets output files\n")
	print_help(parser)
	stop()
} else if (is.null(opt$outFile)) {
	cat("Need output file\n")
	print_help(parser)
	stop()
} else {
	facetsFiles <- arguments$args
}


Df <- data.frame()
for (facetsFile in facetsFiles) {
    load(facetsFile)
    tumorName <- facetsFile %>% sub('.*/', '', .) %>% sub('_.*', '', .) %>% sub('\\..*', '', .)
    normalName <- facetsFile %>% sub('.*/', '', .) %>% sub('^.*_', '', .) %>% sub('\\..*', '', .)
    n <- paste(tumorName, normalName, sep = '_')
    Df[n, 'tumorName'] <- tumorName
    if (tumorName != normalName) {
        Df[n, 'normalName'] <- normalName
    }
    Df[n, 'purity'] <- fit$purity
    Df[n, 'ploidy'] <- fit$ploidy
    Df[n, 'dipLogR'] <- fit$dipLogR
}
Df <- mutate(Df, bad = purity <= 0.3 | is.na(purity))

write.table(Df, file = opt$outFile, sep = '\t', quote = F, row.names = F)
