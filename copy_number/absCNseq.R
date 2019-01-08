#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(absCNseq))
suppressPackageStartupMessages(library(VariantAnnotation))

optionList <- list(
	make_option(c('-t', '--seqType'), action='store', default = "WES", help = 'sequence type (WES or WGS) [Default %default]'),
	make_option(c('-n', '--tumorName'), action='store', default = NULL, help = 'name of the tumor sample. By default, derive from file name'),
	make_option('--genome', action='store', default = 'b37', help = 'genome [default %default]'),
	make_option(c('-o', '--outPrefix'), action='store', default = NULL, help = 'output prefix'))
posArgs <- c('varscanSegFile', 'snvFile')
parser <- OptionParser(usage = paste('%prog [options]', paste(posArgs, collapse=' ')),  option_list=optionList)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args) != length(posArgs)) {
	print_help(parser)
	print(arguments$args)
	stop('Incorrect number of required positional arguments')
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} else {
	cmdArgs <- arguments$args
	for (i in 1:length(cmdArgs)){
		assign(posArgs[i], cmdArgs[i])
	}
}


segData <- read.table(varscanSegFile, sep = '\t', header = T, stringsAsFactors = F)
segRle <- rle(segData$Segmented)
segData <- transform(segData, segId = as.factor(rep(1:length(segRle$values), segRle$lengths)))
segData <- transform(segData, length = End - Start)

chrom <- tapply(segData$Chrom, segData$segId, function(x) x[1])
start <- tapply(segData$Start, segData$segId, function(x) x[1])
end <- tapply(segData$End, segData$segId, function(x) x[length(x)])
effSegLen <- tapply(segData$length, segData$segId, sum)
normRatio <- tapply(segData$Segmented, segData$segId, function(x) x[1])

absSegData <- data.frame(chrom = chrom, loc.start = start, loc.end = end, eff.seg.len = effSegLen, normalized.ratio = normRatio)

if (is.null(opt$tumorName)) {
    tumor <- sub('.*/', '', sub('_.*', '', snvFile))
} else {
    tumor <- opt$tumorName
}
if (grepl('\\.txt$', snvFile, perl = T)) {
    snvData <- read.table(snvFile, sep = '\t', header = T, comment.char = '', stringsAsFactors = F)
    absSnvData <- with(snvData, data.frame(chrom = X.CHROM, position = POS, tumor_var_freq = snvData[,paste(tumor, ".FA", sep = "")]))
} else if (grepl('\\.vcf$', snvFile, perl = T)) {
    vcf <- readVcf(snvFile, opt$genome)
    vaf <- sapply(geno(vcf)$AD[,tumor], function (x) x[2] / sum(x))
    absSnvData <- data.frame(chrom = as.vector(seqnames(vcf)), position = start(rowRanges(vcf)), tumor_var_freq = vaf)
    absSnvData <- subset(absSnvData, !is.na(tumor_var_freq))
} else {
    cat("snv format unreadable\n")
    q(save = 'no', status = 1)
}


if (opt$seqType == "WES") {
    min.seg.len <- 100
} else if (opt$seqType == "WGS") {
    min.seg.len <- 3000
}
res <- grid.search.alpha(absSegData, absSnvData, min.seg.len = min.seg.len)

fn <- paste(opt$outPrefix, ".absCN.txt", sep = "")
write.table(res$absCN, file = fn, sep = '\t', quote = F)

fn <- paste(opt$outPrefix, ".absSNV.txt", sep = "")
write.table(res$absSNV, file = fn, sep = '\t', quote = F)
