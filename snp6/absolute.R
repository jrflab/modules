#!/usr/bin/env Rscript
# Run absolute

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library(ABSOLUTE));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--disease", default = 'breastcancer', help = "disease [default %default]"),
        make_option("--platform", default = "SNP_6.0", help = "platform [default %default]"),
        make_option("--tumour", default = NULL, help = "tumour sample name"),
        make_option("--mafFile", default = NULL, help = "MAF file"),
        make_option("--minMutAF", default = NULL, help = "Minimum Mutation Allele Frequency"),
        make_option("--resultsDir", default = NULL, help = "results directory"),
        make_option("--outPrefix", default = NULL, help = "output prefix")
        )
parser <- OptionParser(usage = "%prog segDat.Rdata", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$resultsDir)) {
    cat("Need results dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$tumour)) {
    cat("Need tumour sample name\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need hapseg data file\n");
    print_help(parser);
    stop();
}


fn <- arguments$args[1];
RunAbsolute(seg.dat.fn = fn, output.fn.base = opt$outPrefix,
    sigma.p=0, max.sigma.h=0.02,
    min.ploidy=0.95, max.ploidy=10, primary.disease=opt$disease,
    platform=opt$platform, sample.name=opt$tumour,
    results.dir=opt$resultsDir,
    maf.fn = opt$mafFile, min.mut.af=opt$minMutAF,
    max.as.seg.count=1500, copy_num_type="allelic",
    max.neg.genome=0, max.non.clonal=0,
    verbose=TRUE)
