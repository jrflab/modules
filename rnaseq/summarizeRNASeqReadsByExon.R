#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

optionList <- list(
	make_option('--genome', action='store', default = 'b37', help = 'genome to use [%default]'),
	make_option(c('-o', '--outFile'), action='store', default = NULL, help = 'output file'))
posArgs <- c('bamFile')
parser <- OptionParser(usage = paste('%prog [options]', paste(posArgs, collapse=' ')),  option_list=optionList)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args) != length(posArgs)) {
	print_help(parser)
	print(arguments$args)
	stop('Incorrect number of required positional arguments')
} else if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else {
	cmdArgs <- arguments$args
	for (i in 1:length(cmdArgs)){
		assign(posArgs[i], cmdArgs[i])
	}
    outFile <- opt$outFile
}

#For debugging
if (FALSE){
	opt <- list('addChr' = TRUE, 'geneListFile' = NULL)
	txdbFile <- '~/ensg69.biomart.13012013.sqlite'
	bamFile <- '~/share/data/DLBCL/WTSS/bam/HS0653.bam'
	outFile <- 'tmp.txt'
}
#txdb <- makeTranscriptDbFromBiomart( biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
#saveFeatures(txdb, '~/share/reference/ensg69.biomart.2014-02-21.sqlite')
#txdb <- makeTranscriptDbFromUCSC(genome = 'hg19', tablename = 'ensGene')

print("Loading txdb ")
if (opt$genome == "b37" || opt$genome == "hg19" || opt$genome == "GRCh37") {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
} else {
    cat("Unsupported genome\n")
    print_help(parser);
    stop();
}

print('... Done')

allExons <- exons(txdb, columns = c('gene_id', 'exon_id', 'exon_name'))

print('Removing chr from chromosome names')
newSeqNames <- sub('chr', '', seqlevels(allExons))
names(newSeqNames) <- seqlevels(allExons)
allExons <- renameSeqlevels( allExons, newSeqNames )

cat("Reading", bamFile, " ... ")
si <- seqinfo(BamFile(bamFile));
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100));
scf <- scanBamFlag( isDuplicate = FALSE ) # remove duplicate reads
reads <- readGappedReads( bamFile, param = ScanBamParam( which = gr, flag = scf ) ); # grab reads in specific region cat("Finished\n")
#reads <- GRanges(seqnames = rname(reads), ranges = IRanges(start = start(reads), end = end(reads)), strand = rep('*', length(reads)));
cat('Finished\n')

print('Count raw exon read counts ...')
#countsForExons <- countOverlaps(allExons, reads);
summarizedExpt <- summarizeOverlaps(allExons, reads)
countsForExons <- as.numeric( assays(summarizedExpt)$counts )
names(countsForExons) <- rownames(summarizedExpt)
print('... Done')

print('Generating expression values ...')
numBases <- width(allExons)
numKBases <- numBases / 1000
millionsMapped <- sum(countsForExons) / 10^6
rpm <- countsForExons / millionsMapped
rpkm <- rpm / numKBases
print('... Done')

print('Retrieving annotation data ...')
annotDf <- values(allExons)
print('...Done')

exonsReadDf <- data.frame( geneID = sapply(annotDf[, 'gene_id'], '[[', 1), exonID = annotDf[, 'exon_id'], exonName = annotDf[, 'exon_name'], exonCount = countsForExons, exonRPM = rpm, exonRPKM = rpkm, stringsAsFactors = FALSE )

print(paste('Writing data to', outFile))
write.table(exonsReadDf, file = outFile, sep = '\t', quote = F, row.names=F)
print('...Done')
