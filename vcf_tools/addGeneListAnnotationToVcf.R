#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("VariantAnnotation"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
}

optList <- list(
        make_option("--genome", default = 'b37', help = "genome build [default %default]"),
        make_option("--geneBed", default = NULL, type = "character", action = "store", help ="input gene beds, comma separated (required)"),
        make_option("--name", default = 'genelist', type = "character", action = "store", help ="annotation names, comma separated (default %default)"),
        make_option("--outFile", default = NULL, type = "character", action = "store", help ="output file (required)"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input vcf file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$geneBed)) {
    cat("Need input gene bed\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file name\n\n")
    print_help(parser);
    stop();
}

Names <- unlist(strsplit(opt$name, ','))
beds <- unlist(strsplit(opt$geneBed, ','))
names(beds) <- Names
if (length(Names) != length(beds)) {
    cat("Number of names != number of bed files")
    print_help(parser)
    stop()
}

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

cat('Reading target bed(s) ... ')
geneLists <- list()
for (n in names(beds)) {
    geneLists[[n]] <- import(beds[[n]])
}
cat('done\n')

vcf <- readVcf(fn, genome = opt$genome)
if (nrow(vcf) > 0) {
    # replace header
    for (n in names(beds)) {
        newInfo <- DataFrame(Number = 0, Type = "Flag", Description = paste(n, ": variant is in gene list", sep = ''), row.names = n)
        info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
        ol <- findOverlaps(rowRanges(vcf), geneLists[[n]], select = 'first')
        info(vcf)[,n] <- !is.na(ol)
    }

    cat("Writing to", opt$outFile, "... ")
    writeVcf(vcf, out)
    cat("done\n")
} else {
    cat("No entries, creating empty vcf file\n")
    vcf <- readVcf(fn, genome = opt$genome)
    writeVcf(vcf, out)
    # fix the empty contig lines
    cmd <- paste("sed -i '/^##contig/d'", outfn)
    system(cmd)
}

close(out)
