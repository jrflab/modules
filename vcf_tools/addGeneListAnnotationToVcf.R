#!/usr/bin/env Rscript
# add gene list annotation column to vcf file
# usage: Rscript addGeneListAnnotationToVcf.R --geneBed [gene lists] --name [gene list names] [vcf file]
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
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

Names <- strsplit(opt$name, ',')
beds <- strsplit(opt$geneBed, ',')
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

cat('Indexing vcf ... ')
temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
cat('done\n')

tab <- TabixFile(zipped, idx, yieldSize = 2000)
open(tab)

cat('Processing vcf by chunk\n')
i <- 1
while(nrow(vcf <- readVcf(tab, genome = opt$genome))) {
    # replace header
    for (n in names(beds)) {
        newInfo <- DataFrame(Number = 0, Type = "Flag", Description = paste(n, ": variant is in gene list", sep = ''), row.names = n)
        info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
        ol <- findOverlaps(rowData(vcf), geneLists[[n]], select = 'first')
        info(vcf)[,n] <- !is.na(ol)
    }

    cat(paste('Chunk', i, "\n"))
    cat("Appending vcf chunk to", opt$outFile, "... ")
    writeVcf(vcf, out)
    cat("done\n")
    i <- i + 1
}

if (i == 1) {
    cat("No entries, creating empty vcf file\n")
    vcf <- readVcf(fn, genome = opt$genome)
    writeVcf(vcf, out)
}
close(tab)
close(out)
