#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("TEQC"))

optList <- list(
                make_option("--ref", default = "hg19", help ="Reference genome [default %default]"),
                make_option("--offset", default = 0, help ="target offset [default %default]"),
                make_option("--outFile", default = NULL, type = "character", action = "store", help ="Output file (required)"))

parser <- OptionParser(usage = "%prog [options] bamFile targetsFile", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need reads file and targets file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

print("reading bam file...")
reads <- get.reads(files[1], filetype = "bam");
x <- sapply(reads, function(x) nrow(x) == 0)
if (any(x)) {
    reads <- reads[-which(x)];
}

print("collapsing reads to pairs...")
readpairs <- reads2pairs(reads)

print("reading targets file...")
targets <- get.targets(files[2]);
baits <- get.baits(files[2]);

print("calculating fraction on targets...")
ft <- fraction.target(targets, genome = opt$ref)
fr <- fraction.reads.target(readpairs, targets, Offset = opt$offset)
enrichment <- fr/ft

print("calculating coverage...")
Coverage <- coverage.target(reads, targets, Offset = opt$offset)
avgcov <- data.frame(round(Coverage$avgTargetCoverage, 2), 
                     round(Coverage$targetCoverageSD, 2), matrix(Coverage$targetCoverageQuantiles, 
                                                                 ncol = 5))
names(avgcov) <- c("avgTargetCoverage", "targetCoverageSD", 
                   paste(names(Coverage$targetCoverageQuantiles), "quantile"))

print("counting reads per target...")
targetcov0 <- Coverage$targetCoverages
targetcov <- readsPerTarget(reads, targetcov0, Offset = opt$offset)

print("counting reads per bait...")
covercounts.baits <- RleList()
baitcov <- NULL
for (chr in names(baits)) {
    cov.chr <- Coverage$coverageAll[[chr]]
    ir.chr <- ranges(baits)[[chr]]
    tmp <- lapply(ir.chr, function(x) seqselect(cov.chr, x))
    avgcov <- sapply(tmp, mean)
    baitcov <- c(baitcov, avgcov)
    cov.chr <- seqselect(cov.chr, reduce(ir.chr))
    covercounts.baits <- c(covercounts.baits, RleList(cov.chr))
}
avgcov <- mean(as.integer(unlist(covercounts.baits)))
baitcov.norm <- baitcov/avgcov

save(targetcov, baitcov, baitcov.norm, avgcov, Coverage, fr, ft, enrichment, targets, baits, file = opt$outFile)
