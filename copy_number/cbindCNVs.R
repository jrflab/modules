#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));
suppressPackageStartupMessages(library("optparse"));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--ensemblTxdb", default = NULL, help = "Ensembl TxDb SQLite"),
        make_option("--outDir", default = NULL, help = "output dir")
        )
parser <- OptionParser(usage = "%prog [CNV files]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outDir)) {
    cat("Need output directory\n");
    print_help(parser);
    stop();
} else if (is.null(opt$ensemblTxdb)) {
    cat("Need ensembl txdb sqlite file\n");
    print_help(parser);
    stop();
}

#arguments <- paste(samples, '.bam_CNVs', sep = '')
# files <- list.files(pattern = "*_CNVs")
files <- arguments$args

grs <- list()
tables <- list()
samples <- c()
for (arg in files) {
    s <- sub('\\..*', '', arg)
    s <- sub('.*/', '', s)
    samples <- c(samples, s)
    d <- read.table(file = arg, sep = '\t', header = F, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start", "end")
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, end = posns$end))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
x <- as.vector(seqnames(gr))
x[x == "X"] <- 23
x[x == "Y"] <- 24
x[x == "MT"] <- 25
x <- as.integer(x)
oo <- order(x, start(gr))
gr <- gr[oo, ]
for (s in samples) {
    mcols(gr)[, s] <- rep(2, length(gr))
}

for (s in samples) {
    overlaps <- findOverlaps(grs[[s]], gr)
    mcols(gr[subjectHits(overlaps), ])[[s]] <- mcols(grs[[s]][queryHits(overlaps), ])$copynum
}


X <- as.matrix(mcols(gr))
X[X > 3] <- 3
#x <- round(width(gr) / 50000)
#Z <- apply(X, 2, rep, times = x)
#sn <- rep(as.vector(seqnames(gr)), times = x)
#chrstart <- seqnames(gr)

#txdb <- makeTranscriptDbFromUCSC('hg19', tablename = 'knownGene')

# annotate positions with genes
#txdb <- makeTranscriptDbFromBiomart('ensembl', 'hsapiens_gene_ensembl')
txdb <- loadDb(opt$ensemblTxdb)
txs <- transcriptsBy(txdb, by = "gene")
overlaps <- findOverlaps(gr, txs)
ensids <- names(txs[subjectHits(overlaps)])
x <- tapply(ensids, queryHits(overlaps), function(x) {
    egids <- sapply(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]);
    if (sum(!is.na(egids) > 0)) {
        genes <- sapply(mget(egids[!is.na(egids)], org.Hs.egSYMBOL), function(x) x[1]);
        paste(genes, collapse = "|");
    } else {
        ""
    }
})
mcols(gr)$Genes <- rep("", length(gr))
mcols(gr)[unique(queryHits(overlaps)), 'Genes'] <- x

# recurrent CNVs by genomic ranges
# positions not diploid in more than one sample
d <- as.data.frame(gr)
x <- rowSums(d[,samples] != 2) > 1
fn <- paste(opt$outDir, '/recurrent_cnv.txt', sep = '')
write.table(d[x,], file = fn, sep = "\t", row.names = F, quote = F)

genes <- unique(unlist(sapply(mcols(gr)$Genes, strsplit, '\\|')))
geneCN <- matrix(2, nrow = length(genes), ncol = length(samples), dimnames = list(Gene = genes, Sample = samples))
for (i in 1:length(gr)) {
    gs <- unlist(strsplit(mcols(gr)[i, "Genes"], '\\|'))
    x <- as.integer(as.data.frame(mcols(gr))[i, samples])
    if (sum(x != 2) > 1) {
        for (g in gs) {
            geneCN[g, samples][x != 2] <- x[x != 2]
        }
    }
}

# recurrent CNVs by gene
# genes not diploid in more than one sample
x <- rowSums(geneCN != 2) > 1
geneRecurrentCNV <- geneCN[x, ]
fn <- paste(opt$outDir, '/gene_recurrent_cnv.txt', sep = '')
write.table(geneRecurrentCNV, file = fn, sep = "\t", quote = F)

x <- rowSums(geneCN > 2) > 1
geneRecurrentGainCNV <- geneCN[x, ]
fn <- paste(opt$outDir, '/gene_recurrent_cnv_gain.txt', sep = '')
write.table(geneRecurrentGainCNV, file = fn, sep = "\t", quote = F)

x <- rowSums(geneCN < 2) > 1
geneRecurrentLossCNV <- geneCN[x, ]
fn <- paste(opt$outDir, '/gene_recurrent_cnv_loss.txt', sep = '')
write.table(geneRecurrentLossCNV, file = fn, sep = "\t", quote = F)

