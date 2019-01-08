#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ExomeCNV"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"),
        make_option("--tumor", default = NULL, type = "character", action = "store", help ="tumor BAF file (required)"),
        make_option("--normal", default = NULL, type = "character", action = "store", help ="normal BAF file (not required)"),
        make_option("--lohMethod", default = "two.sample.fisher", type = "character", action = "store", help ="LoH calling method (default: %default)"),
        make_option("--cbsLohMethod", default = "two.sample.fisher", type = "character", action = "store", help ="CBS LoH calling method (default: %default)"),
        make_option("--alpha", default = 1e-6, action = "store", help ="LOH alpha (default: %default)"))

parser <- OptionParser(usage = "%prog [options] ", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$tumor)) {
    cat("Need tumor BAF file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
}

tumorName <- sub('_.*', '', sub('\\..*', '', sub('.*/', '', opt$tumor)))

if (!is.null(opt$normal)) {
    normalName <- sub('.*_', '', sub('\\..*', '', sub('.*/', '', opt$normal)))
        cat("Reading", normalName, ":", opt$normal, "\n")
        normal <- read.delim(opt$normal, header = T)
} else {
    normal <- NULL
        opt$lohMethod <- "only.tumor"
        opt$cbsLohMethod <- "only.tumor"
}

    cat("Reading", tumorName, ":", opt$tumor, "\n")
tumor <- read.delim(opt$tumor, header = T)


    if (!any(grepl('chr', normal$chr))) {
        if (!is.null(opt$normal)) {
            normal$chr <- sub('^', 'chr', normal$chr)
        }
        tumor$chr <- sub('^', 'chr', tumor$chr)
    }


if (!is.null(opt$normal)) {
    x <- apply(normal, 1, function(x) any(is.na(x))) | apply(tumor, 1, function(x) any(is.na(x)))
        normal <- normal[!x, ]
}
tumor <- tumor[!x, ]


cat("Analyzing LOH using", opt$lohMethod, "\n")
    cat("alpha:", opt$alpha, "\n")
eLOH = LOH.analyze(normal = normal, tumor = tumor, alpha = opt$alpha, method = opt$lohMethod)

    cat("Merging segments\n")
loh = multi.LOH.analyze(normal = normal, tumor = tumor, all.loh.ls = list(eLOH), test.alpha = opt$alpha, method = opt$cbsLohMethod, sdundo = c(1.5,1.5), alpha = c(0.00001,0.00001))
    prefix <- paste(opt$outPrefix, sep = "")
    cat("Writing output (prefix: ", opt$outPrefix, ")\n", sep = "")
    fn <- paste(opt$outPrefix, '.loh.txt', sep = '')
    write.table(loh, file = fn, quote = F, sep = '\t', row.names = F, col.names = T)

    fn <- paste(opt$outPrefix, ".loh.png", sep = "")
    cat("Plotting to", fn, "\n")
    png(filename = fn, res = 70, width = 2000,
            height = 1200, pointsize = 16, type = "cairo-png")
    do.plot.loh(loh, normal = normal, tumor = tumor, method = opt$cbsLohMethod, plot.style = "dev")
dev.off()
