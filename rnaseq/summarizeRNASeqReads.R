#!/usr/bin/env Rscript
# Description: This script is used to generate sum reads over 1) the entire gene, 2) only the exons on a gene, 3) only the introns of a gene. It will also calculate a gene RPM/RPKM based on the exon and introns.
# Description: This script is used to generate sum reads over a gene returning the raw read count across the whole gene and also across only the exons in the gene. Also RPKM values are generated.
# Authors: Raymond Lim and Fong Chun Chan <fongchunchan@gmail.com>
suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(rbamtools));
suppressPackageStartupMessages(library(Rsamtools));
suppressPackageStartupMessages(library(GenomicAlignments));
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene));
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))

optionList <- list(
	#make_option(c('-a', '--addChr'), action='store_true', default = TRUE, help = 'Set the flag to add chr as a prefix to each seqlevel. Necessarily when the bamFile has read containing the chr prefix and the txdbFile does not [%default]'),
	make_option(c('-i', '--intronListFile'), action='store', default = NULL, help = 'Set a file containing intronIDs to include in the summarization [%default]'),
	make_option(c('-g', '--genome'), action='store', default = 'b37', help = 'genome to use [%default]'),
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


#For debugging
if (FALSE){
#	opt <- list( 'addChr' = TRUE, 'intronListFile' = '~/dlbcl_snp6/paper.analysis.code/references/introns/introns.non.overlapping.with.exons.txt', 'intronWindow' = 50 )
	opt <- list( 'addChr' = TRUE, 'countMethod' = 'summarizeOverlaps', 'intronListFile' = NULL, 'intronWindow' = NULL )
	txdbFile <- '~/ensg69.biomart.13012013.sqlite'
#	txdbFile <- '~/hg19_ensGene.06022012.sqlite'
	bamFile <- '~/share/data/DLBCL/WTSS/tophat/bam/HS0751.bam'
	outFile <- 'tmp.txt'
}

#txdb <- makeTranscriptDbFromBiomart( biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl' )
#saveDb(txdb, '~/GRCm38.08032013.sqlite')
#txdb <- makeTranscriptDbFromUCSC(genome = 'hg19', tablename = 'ensGene')

# Description: This function is used to get the raw reads over genomic features
# Inputs: 
#	1) features : A GRanges object
#	2) reads: A GappedAlignments object generated from readBamGappedAlignments
#	3) countMethod: Two alternative counting methods can be specified. i) countOverlaps allows for reads to be counted by than once, ii) summarizeOverlaps counts reads only once 
# Outputs:
getCounts <- function( features, reads ){
	summarizedExpt <- summarizeOverlaps( features, reads )
	counts <- as.numeric( assays(summarizedExpt)$counts )
	names(counts) <- rownames(summarizedExpt)
	return(counts)
}

# Description: This function is used to get RPKM values over genomic features
# Inputs: 
#	1) features : A GRanges object
#	2) featureCounts: A numeric vector containing the number of raw counts aligning to the feature
# Outputs:
getExprs <- function( features, featureCounts, feature = 'gene' ){
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
} else if (opt$genome == "mm10" || opt$genome == "GRCm38") {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
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
	intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '') # intronID should be the whole intron genomic coordinate 
	intronList <- scan( opt$intronListFile, what = 'character' )
	intronFlag <- intronIDs %in% intronList
	introns <- introns[ intronFlag ]
	intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '') # intronID should be the whole intron genomic coordinate 
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
si <- seqinfo(BamFile(bamFile));
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100));
scf <- scanBamFlag( isDuplicate = FALSE ) # remove duplicate reads
reads <- readGappedReads( bamFile, param = ScanBamParam( which = gr, flag = scf )); # grab reads in specific region
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
