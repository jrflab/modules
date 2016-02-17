#!/usr/bin/env Rscript
# add somatic sniper AD field

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--genome", default = 'b37', help = "genome build [default %default]"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]")
        )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

cat('Reading vcf header ... ')
# create new header
vcfHeader <- scanVcfHeader(fn)
gen <- apply(as.data.frame(geno(vcfHeader)), 2, as.character)
rownames(gen) <- rownames(geno(vcfHeader))
newGeno <- DataFrame(rbind(gen, AD = c(".", "Integer", "Allelic depth")))
hlist <- header(vcfHeader)
hlist$FORMAT <- newGeno
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)
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
    oldwd <- getwd()
    # replace header
    metadata(vcf)$header <- newVcfHeader
    # pre-populate new info fields with NAs
    newGeno <- geno(vcf)
    newGeno$AD <- aperm(apply(geno(vcf)$DP4, c(1,2), function(x) c(x[1] + x[2], x[3] + x[4])), c(2, 3, 1))
    geno(vcf) <- newGeno

    cat(paste('Chunk', i, "\n"))
    i <- i + 1

    #fix sample genotype order
    x <- which(names(geno(vcf)) == "GT")
    ord <- c(x, (1:length(geno(vcf)))[-x])
    geno(vcf) <- geno(vcf)[ord]

    cat("Appending vcf chunk to", opt$outFile, "... ")
    setwd(oldwd)
    writeVcf(vcf, out)
    cat("done\n")
}

close(tab)
close(out)

