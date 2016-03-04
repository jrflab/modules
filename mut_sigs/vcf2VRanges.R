#!/usr/bin/env Rscript
# convert a vcf file to a VRanges Rdata object
library("optparse");
library("VariantAnnotation");
library("BSgenome.Hsapiens.UCSC.hg19");
library("BSgenome.Mmusculus.UCSC.mm10");
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("reshape")
library("boot")
library("plyr")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("reshape2")
library("SomaticSignatures")
library("foreach")

optList <- list(
                make_option("--genome", default = 'b37', help = "reference genome"),
                make_option("--ignoreFilter", default = F, action = 'store_true', help = "ignore the filter column for vcf files"),
                make_option("--outFile", default = NULL, type = "character", action = "store", help = "output directory"))

parser <- OptionParser(usage = "%prog [options] [vcf file(s)]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

vcfFile <- arguments$args[1]
outFile <- opt$outFile
if (opt$genome == "b37" || opt$genome == "hg19") {
    genome <- BSgenome.Hsapiens.UCSC.hg19
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    genomeName <- 'hg19'
    chromosomes <- c(1:22, "X", "Y")
} else if (opt$genome == "mm10" || opt$genome == "GRCm38") {
    genome <- BSgenome.Mmusculus.UCSC.mm10
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    genomeName <- 'mm10'
    chromosomes <- c(1:19, "X", "Y")
}

txByGenes <- transcriptsBy(txdb, 'gene')

temp <- tempfile()
zipped <- bgzip(vcfFile, temp)
idx <- indexTabix(temp, "vcf")
cat('done\n')

tab <- TabixFile(zipped, idx)
open(tab)

vcf <- readVcf(tab, genomeName)
passIds <- which(rowRanges(vcf)$FILTER == "PASS")
if (nrow(vcf) > 0 && sum(seqnames(vcf) %in% chromosomes) > 0 &&
    sum(isSNV(vcf)) > 0 && (opt$ignoreFilter | length(passIds) > 0)) {
    vcf <- vcf[isSNV(vcf) & seqnames(vcf) %in% chromosomes]
    if (!opt$ignoreFilter) {
        vcf <- vcf[passIds, ]
    }
    s <- sub('\\..*', '', vcfFile)
    s <- sub('.*/', '', s)
    vr <- VRanges(seqnames = seqnames(vcf),
            ranges = ranges(vcf),
            ref = as.character(ref(vcf)),
            alt = sapply(alt(vcf), function(x) as.character(x[1])),
            sampleNames = s)
    seqlevels(vr) <- sub('^M$', 'MT', seqlevels(vr))
    vr <- ucsc(vr)
    vr <- mutationContext(vr, genome, unify = T)
    vr$refalt <- paste(ref(vr), alt(vr), sep = '')

    # query transcript ids
    ol <- findOverlaps(vr, txByGenes)
    subjectStrands <- sapply(txByGenes[subjectHits(ol)], function(x) paste(unique(as.character(strand(x))), collapse = ','))
    queryStrands <- tapply(subjectStrands, queryHits(ol), function(x) paste(unique(x), collapse = ","))
    vr$txStrand <- NA
    vr$txStrand[as.integer(names(queryStrands))] <- queryStrands
    vr$transcribed <- F
    vr$transcribed[is.na(vr$txStrand)] <- NA
    vr$transcribed[vr$refalt %in% c("GA", "GC", "GT", "AC", "AG", "AT") & grepl('\\+', vr$txStrand)] <- T
    vr$transcribed[vr$refalt %in% c("CA", "CG", "CT", "TA", "TC", "TG") & grepl('-', vr$txStrand)] <- T
    save(vr, file = opt$outFile)
} else {
    vr <- NULL
    save(vr, file = opt$outFile)
}
