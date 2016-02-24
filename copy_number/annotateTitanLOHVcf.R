#!/usr/bin/env Rscript
# annotate vcf file with titan seg file

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"));
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GenomicRanges"));

#options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--titanSeg", default = NULL, type = "character", action = "store", help ="targeted titan segment file"),
        make_option("--outFile", default = NULL, type = "character", action = "store", help ="targeted interval bed"),
        make_option("--genome", default = 'b37', type = "character", action = "store", help ="reference genome"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need vcf file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$titanSeg)) {
    cat("Need titan seg file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n\n")
    print_help(parser);
    stop();
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- sub('chr', '', seqlevels(txdb))
tx <- transcripts(txdb)

# read titan file
titanSeg <- read.table(opt$titanSeg, sep = '\t', header = T)
titanGR <- with(titanSeg, GRanges(seqnames = Chromosome, ranges = IRanges(start = Start_Position.bp., end = End_Position.bp.), CN = Copy_Number, minorCN = MinorCN, majorCN = MajorCN, call = TITAN_call, medianRatio = Median_Ratio, medianLogR = Median_logR))
ol <- findOverlaps(tx, titanGR, select = 'first')
tx$CN[!is.na(ol)] <- titanGR$CN[ol[!is.na(ol)]]
tx$minorCN[!is.na(ol)] <- titanGR$minorCN[ol[!is.na(ol)]]
tx$majorCN[!is.na(ol)] <- titanGR$majorCN[ol[!is.na(ol)]]
tx$titanCall[!is.na(ol)] <- titanGR$call[ol[!is.na(ol)]]
tx$medianRatio[!is.na(ol)] <- titanGR$medianRatio[ol[!is.na(ol)]]
tx$medianLogR[!is.na(ol)] <- titanGR$medianLogR[ol[!is.na(ol)]]
titanTx <- transcriptsByOverlaps(txdb, titanGR)
ol <- findOverlaps(titanTx, titanGR, select = 'first')
mcols(titanTx) <- cbind(mcols(titanTx), mcols(titanGR[ol]))

# output file
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

#input file
vcfFile <- arguments$args[1]
temp <- tempfile()
zipped <- bgzip(vcfFile, temp)
idx <- indexTabix(temp, "vcf")
tab <- TabixFile(zipped, idx, yieldSize = 8000)

open(tab)
while(nrow(vcf <- readVcf(tab, opt$genome))) {
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Titan transcript copy number", row.names = "titanCN"))
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Titan transcript minor allele copy number", row.names = "titanMinorCN"))
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer", Description = "Titan transcript major allele copy number", row.names = "titanMajorCN"))
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String", Description = "Titan transcript call", row.names = "titanCall"))
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String", Description = "Titan transcript median ratio", row.names = "titanMedianRatio"))
    info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String", Description = "Titan transcript median logR", row.names = "titanMedianLogR"))

    seqlevels(vcf) <- sub('chr', '', seqlevels(vcf))
    ol <- findOverlaps(rowRanges(vcf), titanTx, select = 'first')
    info(vcf)$titanCN[!is.na(ol)] <- titanTx$CN[ol[!is.na(ol)]]
    info(vcf)$titanMinorCN[!is.na(ol)] <- titanTx$minorCN[ol[!is.na(ol)]]
    info(vcf)$titanMajorCN[!is.na(ol)] <- titanTx$majorCN[ol[!is.na(ol)]]
    info(vcf)$titanCall[!is.na(ol)] <- as.character(titanTx$call[ol[!is.na(ol)]])
    info(vcf)$titanMedianRatio[!is.na(ol)] <- titanTx$medianRatio[ol[!is.na(ol)]]
    info(vcf)$titanMedianLogR[!is.na(ol)] <- titanTx$medianLogR[ol[!is.na(ol)]]
    writeVcf(vcf, out)
}
close(tab)
close(out)
