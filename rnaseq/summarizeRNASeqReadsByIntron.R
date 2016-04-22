#!/usr/bin/env Rscript
# Description: This script is used to generate intronic counts. You can pass in an optional parameter, intronWindow, that restrict the intronic window to a certain number of bases. Also, can pass in an optional intron list file which will restrict the summarization to just these introns
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
library("GenomicFeatures")
library("GenomicAlignments");
library("Rsamtools")
library('optparse')

optionList <- list(
	make_option(c('-a', '--addChr'), action='store_true', default = FALSE, help = 'Set the flag to add chr as a prefix to each seqlevel [%default]'),
	make_option(c('-i', '--intronListFile'), action='store', default = NULL, help = 'Set a file containing intronIDs to include in the summarization [%default]'),
	make_option(c('-d', '--txdb'), action='store', default = NULL, help = 'ensembl transcript database'),
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
} else if (is.null(opt$txdb)) {
    cat("Need ensembl transcript database\n");
    print_help(parser);
    stop();
} else {
	cmdArgs <- arguments$args
	for (i in 1:length(cmdArgs)){
		assign(posArgs[i], cmdArgs[i])
	}
    txdbFile <- opt$txdb
    outFile <- opt$outFile
}

#For debugging
if (FALSE){
	opt <- list('addChr' = TRUE, 'intronListFile' = '~/dlbcl_snp6/paper.analysis.code/references/introns/introns.non.overlapping.with.exons.txt', 'intronWindow' = 50 )
	txdbFile <- '~/ensg69.biomart.13012013.sqlite'
#	txdbFile <- '~/hg19_ensGene.06022012.sqlite'
	bamFile <- '~/share/data/DLBCL/WTSS/bam/HS0653.bam'
	outFile <- 'tmp.txt'
}

#txdb <- makeTranscriptDbFromBiomart( biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
#saveFeatures(txdb, '~/ensg69.biomart.13012013.sqlite')
#txdb <- makeTranscriptDbFromUCSC(genome = 'hg19', tablename = 'ensGene')

cat("Loading", txdbFile, " ... ")
txdb <- loadDb(txdbFile)
print('... Done')

print('Getting introns by genes ...')
intronsByTx <- intronsByTranscript(txdb, use.names = TRUE )
introns <- unlist(intronsByTx)
introns <- introns[ !duplicated(introns) ]

if ( opt$addChr ){
	print('Adding the chr prefix ...')
	newSeqNames <- paste('chr', seqlevels(introns), sep = '')
	names(newSeqNames) <- seqlevels(introns)
	introns <- renameSeqlevels( introns, newSeqNames )
	print('... Done')
} 
intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '') # intronID should be the whole intron genomic coordinate 

if ( !is.null(opt$intronListFile) ){
	print('Filtering the intron list ...')
	intronList <- scan( opt$intronListFile, what = 'character' )
	intronIDs <- gsub( 'chrX', '23', intronIDs)
	intronIDs <- gsub( 'chrY', '24', intronIDs)
	intronFlag <- intronIDs %in% intronList
	introns <- introns[ intronFlag ]
	intronIDs <- paste(seqnames(introns), ':', start(introns), '-', end(introns), sep = '') # intronID should be the whole intron genomic coordinate 
	print('... Done')
}

if ( !is.null(opt$intronWindow) ){
	print(paste('Restricting the intronic window from intron start +', opt$intronWindow))
	end(introns) <- start(introns) + opt$intronWindow
	print('... Done')
}
print('... Done')

cat("Reading", bamFile, " ... ")
si <- seqinfo(BamFile(bamFile));
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100));
scf <- scanBamFlag( isDuplicate = FALSE ) # remove duplicate reads
reads <- readGappedReads( bamFile, param = ScanBamParam( which = gr, flag = scf ) ); # grab reads in specific region
#reads <- GRanges(seqnames = rname(reads), ranges = IRanges(start = start(reads), end = end(reads)), strand = rep('*', length(reads)));
cat('Finished\n')

print('Count raw intron read counts ...')
#countsForIntrons <- countOverlaps(introns, reads);
summarizedExpt <- summarizeOverlaps(introns, reads)
countsForIntrons <- as.numeric( assays(summarizedExpt)$counts )
names(countsForIntrons) <- rownames(summarizedExpt)
print('... Done')

print('Generating expression values ...')
numBases <- width(introns)
numKBases <- numBases / 1000
millionsMapped <- sum(countsForIntrons) / 10^6
rpm <- countsForIntrons / millionsMapped
rpkm <- rpm / numKBases
print('... Done')

print( paste('Saving results to', outFile) )
map <- select(txdb, keys = names(introns), keytype='TXNAME', cols='GENEID')
geneIDs <- map$GENEID
summarizedReads <- data.frame(geneID = geneIDs, intronID = intronIDs, intronName = intronIDs, chr = as.vector(seqnames(introns)), start = start(introns), stop = end(introns), intronCount = countsForIntrons, intronRPM = rpm, intronRPKM = rpkm)
summarizedReads[, 'chr'] <- gsub('X', 23, summarizedReads[, 'chr'])
summarizedReads[, 'chr'] <- gsub('Y', 24, summarizedReads[, 'chr'])
summarizedReads <- summarizedReads[order(summarizedReads$chr, summarizedReads$start, summarizedReads$stop), ] #Order by the gene chr, gene start and then gene end
summarizedReads <- summarizedReads[, c('geneID', 'intronID', 'intronName', 'intronCount', 'intronRPM', 'intronRPKM')]
write.table(summarizedReads, file = outFile, sep = '\t', quote = F, row.names=F)
print('... Done')
