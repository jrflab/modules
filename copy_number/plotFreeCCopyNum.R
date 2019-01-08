#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outPrefix", default = NULL, help = "Output file prefix"),
                make_option("--includeChrY", action = 'store_true', default = F, help = "include Y-chromosome"),
                make_option("--centromereTable", help = "Centromere position table"));

parser <- OptionParser(usage = "%prog [options] [list of ratio.txt files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input controlFreeC ratio.txt files\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}


grs <- list()
tables <- list()
samples <- c()
for (f in files) {
    s <- sub('\\..*', '', f)
    s <- sub('.*/', '', s)
    samples <- c(samples, s)
    d <- read.table(file = f, sep = '\t', header = T, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], width = 50000), medianratio = d[,4], copynum = d[,5])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start" )
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, width = 50000))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
if (!opt$includeChrY) {
    gr <- gr[seqnames(gr) != "Y"]
}
x <- as.vector(seqnames(gr))
x[x == "X"] <- 23
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

if (!is.null(opt$centromereTable)) {
    centromerePosns <- subset(read.table(opt$centromereTable, sep = '\t'), V8 == "centromere", as.is = T)
    centromerePosns <- centromerePosns[, 2:4]
    centromerePosns[, 1] <- sub('chr', '', centromerePosns[, 1])
    centromerePosns[centromerePosns[,1] == "X", 1] <- 23
    centromerePosns <- centromerePosns[-which(centromerePosns[,1] == "Y"), ]
    oo <- order(as.integer(centromerePosns[,1]))
    centromerePosns <- centromerePosns[oo, ]
    centromerePosns[centromerePosns[,1] == 23, 1] <- "X"
    centromereGR <- GRanges(seqnames = centromerePosns[, 1], ranges = IRanges(centromerePosns[, 2], end = centromerePosns[, 3]))

    cmPos <- nearest(centromereGR, gr)
    cmPos <- (cmPos - 1) / (length(gr) - 1)
}


chrPos <- (start(seqnames(gr)) - 1) / (length(gr) - 1)

X <- as.matrix(mcols(gr))
X[X == 4] <- 3
X[X > 4] <- 4
#X <- X[,-4]
#x <- round(width(gr) / 50000)
#Z <- apply(X, 2, rep, times = x)
#sn <- rep(as.vector(seqnames(gr)), times = x)
#chrstart <- seqnames(gr)


cols <- c('Red', 'Darkred', 'white', 'Darkgreen', 'green')
rng <- range(X)
cols <- cols[(rng[1]:rng[2])+1]
fn <- paste(opt$outPrefix, ".large.png", sep = '')
png(fn, height = 1000, width = 30000, type = 'cairo-png')
par(mar = c(5,10,5,5))
image(X, col = cols, axes = F)
legend('top', legend = c("Del", "Loss", "Neutral", "Gain", "Amp"), fill = cols, horiz = T)
axis(1, at = chrPos, labels = as.character(runValue(seqnames(gr))), cex.axis = 1.5)
abline(v = chrPos, col = 'grey')
if (!is.null(opt$centromereTable)) {
    abline(v = cmPos, lty = 2, col = 'grey')
}
if (ncol(X) == 1) {
    axis(2, at = 0, labels = colnames(X), cex.axis = 1.5, las = 2)
} else {
    axis(2, at = 0:(ncol(X) - 1) / (ncol(X) - 1), labels = colnames(X), cex.axis = 1.5, las = 2)
}
box()
#legend('top', legend = c("deletion", "loss", "neutral", "gain"), fill = cols, horiz = T)
null <- dev.off()

fn <- paste(opt$outPrefix, ".png", sep = '')
png(fn, height = 1000, width = 3000, type = 'cairo-png')
par(mar = c(5,10,5,5))
image(X, col = cols, axes = F)
axis(1, at = chrPos, labels = as.character(runValue(seqnames(gr))), cex.axis = 2)
abline(v = chrPos, col = 'grey')
if (!is.null(opt$centromereTable)) {
    abline(v = cmPos, lty = 2, col = 'grey')
}
if (ncol(X) == 1) {
    axis(2, at = 0, labels = colnames(X), cex.axis = 2, las = 2)
} else {
    axis(2, at = 0:(ncol(X) - 1) / (ncol(X) - 1), labels = colnames(X), las = 2, cex.axis = 2)
}
box()
legend('top', legend = c("Del", "Loss", "Neutral", "Gain", "Amp"), fill = cols, horiz = T)
#legend('top', legend = c("deletion", "loss", "neutral", "gain"), fill = cols, horiz = T)
null <- dev.off()
#write.table(X, file="~/X.txt", sep="\t", col.names=NA, quote=F)
#save(gr, file="~/gr.RData")
