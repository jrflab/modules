#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))


optionList <- list(
	make_option(c('-i', '--intronListFile'), action='store', default = NULL, help = 'Set a file containing intronIDs to include in the summarization [%default]'),
	make_option('--genome', action='store', default = 'b37', help = 'genome to use [%default]'),
	make_option(c('-o', '--outFile'), action='store', default = NULL, help = 'output file'),
	make_option(c('-w', '--intronWindow'), action='store', type = 'integer', default = NULL, help = 'Set the intronic window to be first x bases from the start [%default]')
	)
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
}


if (FALSE) {
	opt <- list( 'addChr' = TRUE, 'countMethod' = 'summarizeOverlaps', 'intronListFile' = NULL, 'intronWindow' = NULL )
	txdbFile <- '~/ensg69.biomart.13012013.sqlite'
	bamFile <- '~/share/data/DLBCL/WTSS/tophat/bam/HS0751.bam'
	outFile <- 'tmp.txt'
}

getCounts <- function( features, reads )
{
	summarizedExpt <- summarizeOverlaps( features, reads )
	counts <- as.numeric( assays(summarizedExpt)$counts )
	names(counts) <- rownames(summarizedExpt)
	return(counts)
}

getExprs <- function( features, featureCounts, feature = 'gene' )
{
	numBases <- sum(width(features))
	numKBases <- numBases / 1000
	millionsMapped <- sum(featureCounts) / 10^6
	rpm <- featureCounts / millionsMapped
	rpkm <- rpm / numKBases
	return( list( 'rpm' = rpm, 'rpkm' = rpkm) )
}

print("Loading txdb ")
if (opt$genome == "b37" || opt$genome == "hg19" || opt$genome == "GRCh37") {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
} else {
    cat("Unsupported genome\n")
    print_help(parser);
    stop();
}

print('... Done')

print('Getting transcripts by gene')
txByGene <- transcriptsBy( txdb, 'gene' )

print('Getting exons by genes ...')
exonsByGene <- exonsBy(txdb, 'gene')
print('... Done')

print('Getting introns by genes ...')
intronsByTx <- intronsByTranscript(txdb, use.names = TRUE )
introns <- unlist(intronsByTx)
introns <- introns[ !duplicated(introns) ]

print('Removing chr from chromosome names')
newSeqNames <- sub('chr', '', seqlevels(txByGene))
names(newSeqNames) <- seqlevels(txByGene)
txByGene <- renameSeqlevels( txByGene, newSeqNames )

newSeqNames <- sub('chr', '', seqlevels(exonsByGene))
names(newSeqNames) <- seqlevels(exonsByGene)
exonsByGene <- renameSeqlevels( exonsByGene, newSeqNames )

newSeqNames <- sub('chr', '', seqlevels(introns))
names(newSeqNames) <- seqlevels(introns)
introns <- renameSeqlevels( introns, newSeqNames )
print('... Done')

if ( !is.null(opt$intronWindow) ){
	print(paste('Restricting the intronic window from intron start +', opt$intronWindow))
	end(introns) <- start(introns) + opt$intronWindow
	print('... Done')
}

if ( !is.null(opt$intronListFile) ){
	print('Filtering the intron list ...')
	intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '')
	intronList <- scan( opt$intronListFile, what = 'character' )
	intronFlag <- intronIDs %in% intronList
	introns <- introns[ intronFlag ]
	intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '')
	print('... Done')
}

txnames <- names(introns)
map <- select(txdb, keys=txnames, keytype='TXNAME', columns='GENEID')
idx <- map$GENEID[!is.na(map$GENEID)]
intronsByGene <- split(introns[!is.na(map$GENEID)], idx)
names(intronsByGene) <- unique(idx)

genes <- c( names(txByGene), names(exonsByGene), names(intronsByGene) )
genes <- unique( genes )
print('... Done')

cat("Reading", bamFile, " ... ")
si <- seqinfo(BamFile(bamFile))
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100))
scf <- scanBamFlag( isDuplicate = FALSE )
reads <- readGappedReads( bamFile, param = ScanBamParam( which = gr, flag = scf ))
cat('Finished\n')

print('Summarizing raw reads over the exon and introns ...')
countsForGenes <- getCounts( txByGene, reads )
countsForExons <- getCounts( exonsByGene, reads )
countsForIntrons <- getCounts( intronsByGene, reads )
cat('Finished\n')

print('Generating expression values ...')
print('...for exons ...')
exonExprsVals <- getExprs( exonsByGene, countsForExons )

print('... for introns ...')
intronExprsVals <- getExprs( intronsByGene, countsForIntrons )
print('... Done')

if (opt$genome == "hg19" || opt$genome == "b37" || opt$genome == "GRCh37") {
    geneSymbols <- sapply(mget(genes, org.Hs.egSYMBOL, ifnotfound = NA), function (x) x[1])
} else {
    geneSymbols <- sapply(mget(genes, org.Mm.egSYMBOL, ifnotfound = NA), function (x) x[1])
}

print( paste('Saving results to', opt$outFile) )
summarizedReads <- data.frame(geneID = genes, gene = geneSymbols, countsByGene = countsForGenes[genes], countsByExon = countsForExons[genes], countsByIntron = countsForIntrons[genes], exonRPM = exonExprsVals[['rpm']][genes], exonRPKM = exonExprsVals[['rpkm']][genes], intronRPM = intronExprsVals[['rpm']][genes], intronRPKM = intronExprsVals[['rpkm']][genes])
write.table(summarizedReads, file = opt$outFile, sep = '\t', quote = F, row.names=F)
print('... Done')
