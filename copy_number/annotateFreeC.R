#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, help = "Output dir"),
                make_option("--includeChrY", action = 'store_true', default = F, help = "include Y chromosome"),
                make_option("--knownVariants", default = NULL, help = "known variants file"),
                make_option("--txdb", default = NULL, help = "txdb"))

parser <- OptionParser(usage = "%prog [options] [list of ratio.txt files]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input controlFreeC CNV files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

#files <- c('freec/AdCCPCT_AdCCPCN/AdCCPCT.bam_CNVs', 'freec/AdCC13T_AdCCPCN/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCCPCN/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCCPCN/AdCC7T.bam_CNVs', 'freec/AdCCPC2T_AdCCPC2N/AdCCPC2T.bam_CNVs') #, 'freec/AdCC13T_AdCCPC2N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCCPC2N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCCPC2N/AdCC7T.bam_CNVs', 'freec/AdCCKHT_AdCCKHN/AdCCKHT.bam_CNVs', 'freec/AdCC13T_AdCCKHN/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCCKHN/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCCKHN/AdCC7T.bam_CNVs', 'freec/AdCC9T_AdCC9N/AdCC9T.bam_CNVs', 'freec/AdCC13T_AdCC9N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC9N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC9N/AdCC7T.bam_CNVs', 'freec/AdCC9repT_AdCC9repN/AdCC9repT.bam_CNVs', 'freec/AdCC13T_AdCC9repN/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC9repN/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC9repN/AdCC7T.bam_CNVs', 'freec/AdCC8T_AdCC8N/AdCC8T.bam_CNVs', 'freec/AdCC13T_AdCC8N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC8N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC8N/AdCC7T.bam_CNVs', 'freec/AdCC6T_AdCC6N/AdCC6T.bam_CNVs', 'freec/AdCC13T_AdCC6N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC6N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC6N/AdCC7T.bam_CNVs', 'freec/AdCC5T_AdCC5N/AdCC5T.bam_CNVs', 'freec/AdCC13T_AdCC5N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC5N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC5N/AdCC7T.bam_CNVs', 'freec/AdCC4T_AdCC4N/AdCC4T.bam_CNVs', 'freec/AdCC13T_AdCC4N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC4N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC4N/AdCC7T.bam_CNVs', 'freec/AdCC3T_AdCC3N/AdCC3T.bam_CNVs', 'freec/AdCC13T_AdCC3N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC3N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC3N/AdCC7T.bam_CNVs', 'freec/AdCC32T_AdCC32N/AdCC32T.bam_CNVs', 'freec/AdCC13T_AdCC32N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC32N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC32N/AdCC7T.bam_CNVs', 'freec/AdCC2T_AdCC2N/AdCC2T.bam_CNVs', 'freec/AdCC13T_AdCC2N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC2N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC2N/AdCC7T.bam_CNVs', 'freec/AdCC1T_AdCC1N/AdCC1T.bam_CNVs', 'freec/AdCC13T_AdCC1N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC1N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC1N/AdCC7T.bam_CNVs', 'freec/AdCC12T_AdCC12N/AdCC12T.bam_CNVs', 'freec/AdCC13T_AdCC12N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC12N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC12N/AdCC7T.bam_CNVs', 'freec/AdCC11T_AdCC11N/AdCC11T.bam_CNVs', 'freec/AdCC13T_AdCC11N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC11N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC11N/AdCC7T.bam_CNVs', 'freec/AdCC10T_AdCC10N/AdCC10T.bam_CNVs', 'freec/AdCC13T_AdCC10N/AdCC13T.bam_CNVs', 'freec/AdCCPC3T_AdCC10N/AdCCPC3T.bam_CNVs', 'freec/AdCC7T_AdCC10N/AdCC7T.bam_CNVs')

grs <- list()
tables <- list()
samples <- c()
for (f in files) {
    s <- sub('\\..*', '', f)
    s <- sub('freec/', '', s)
    s <- sub('/.*', '', s)
    samples <- c(samples, s)
    d <- read.table(file = f, sep = '\t', header = F, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start", "end")
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, end = posns$end))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
if (!opt$includeChrY) {
    gr <- gr[seqnames(gr) != "Y"]
}
x <- as.vector(seqnames(gr))
if (any(x == "X")) {
    x[x == "X"] <- 23
}
if (any(x == "Y")) {
    x[x == "Y"] <- 24
}
if (any(x == "MT")) {
    x[x == "MT"] <- 25
}
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


if (!is.null(opt$txdb)) {
    txdb <- loadDb(opt$txdb)
} else {
    txdb <- makeTranscriptDbFromBiomart('ensembl', 'hsapiens_gene_ensembl')
}

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

fn <- paste(opt$outDir, "/gene_copynum_recurrent.txt", sep = "")
x <- rowSums(geneCN != 2) > 1
geneRecurrentCNV <- geneCN[x, ]
write.table(geneRecurrentCNV, file = fn, sep = "\t", quote = F)

fn <- paste(opt$outDir, "/gene_copynum_recurrent_gain.txt", sep = "")
x <- rowSums(geneCN > 2) > 1
geneRecurrentGainCNV <- geneCN[x, ]
write.table(geneRecurrentGainCNV, file = fn, sep = "\t", quote = F)

fn <- paste(opt$outDir, "/gene_copynum_recurrent_loss.txt", sep = "")
x <- rowSums(geneCN < 2) > 1
geneRecurrentLossCNV <- geneCN[x, ]
write.table(geneRecurrentLossCNV, file = fn, sep = "\t", quote = F)

if (!is.null(opt$knownVariants)) {
    dgv <- read.table(opt$knownVariants, header = T, quote = '', comment.char = '', as.is = T, sep = '\t')
        grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
    dgv.gr <- GRanges(seqnames = dgv$chr, ranges = IRanges(dgv$start, end = dgv$end), id = dgv$variantaccession, observedgains = dgv$observedgains, observedlosses = dgv$observedlosses)
    # take the small known variants (< 500 kb)
    small.dgv.gr <- sort(dgv.gr[width(dgv.gr) <= 500000])
    # at least half overlap
    ols <- findOverlaps(gr, small.dgv.gr, minoverlap = 25000L)
    rs <- ranges(ols, ranges(gr), ranges(small.dgv.gr))
    olf <- width(rs) / width(gr)[queryHits(ols)]
    x <- ols[olf > 0.5]
    mcols(gr)$ID <- rep(".", length(gr))
    ids <- mcols(dgv.gr[subjectHits(x)])$id
    ids <- tapply(ids, queryHits(x), paste, collapse = '|')
    small.gr.index <- which(width(gr) <= 500000)
    xx <- intersect(as.integer(names(ids)), small.gr.index)
    mcols(gr)$ID[xx] <- ids[xx]
    #mcols(gr) <- cbind(mcols(gr), mcols(small.dgv.gr[subjectHits(x)]))
}

fn <- paste(opt$outDir, "/annotated_cnv.txt", sep = "")
write.table(as.data.frame(gr), file = fn, sep = "\t", row.names = F, quote = F)
