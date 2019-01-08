#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("beeswarm"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("TEQC"))
suppressPackageStartupMessages(library("hwriter"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomeGraphs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--covThreshold", default = NULL, help = "coverage hist threshold [default %default]"),
                make_option("--outDir", default = NULL, type = "character", action = "store", help ="Output directory (required)"))

parser <- OptionParser(usage = "%prog [options] RdataFile1 RdataFile2 .. RdataFileN", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need R data file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

avgcovs <- list()
Coverages <- list()
enrichments <- list()
frs <- list()
targetcovs <- list()
baitcovs <- list()
baitcov.norms <- list()

### LOADING R DATA ###
cat("Loading R data: ")
for (f in files) {
    n <- sub('\\..*', '', basename(f))
    cat(n, " ")
    load(f)
    avgcovs[[n]] <- avgcov
    Coverages[[n]] <- Coverage
    enrichments[[n]] <- enrichment
    targetcovs[[n]] <- targetcov
    baitcovs[[n]] <- baitcov
    baitcov.norms[[n]] <- baitcov.norm
    frs[[n]] <- fr
}
cat(' done\n')

avgcov <- do.call('rbind', avgcovs)

targetcov.norms <- lapply(targetcovs, function(x) x[["avgCoverage"]] / mean(x[["nReads"]]))

pg <- openPage("index.html", dirname = opt$outDir)


### COVERAGE PLOTS ###
cat("Generating coverage histogram ... ")
gfn <- paste(opt$outDir, '/coverageHist.png', sep = '')
png(gfn, height = 600, width = 600, type = 'cairo-png')
hist(avgcov, xlab = 'Average Coverage', main = "Average Bait Coverage Distribution")
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)
cat(" done\n")


cat("Generating coverage uniformity plots ... ")
gfn <- paste(opt$outDir, '/coverageUniformity.png', sep = '')
png(gfn, height = 1200, width = 1200, type = 'cairo-png')
n <- ceiling(length(Coverages) / 8)
ltys <- as.integer(gl(n, 8))
cols <- rep(1:8, n)
coverage.uniformity(Coverages[[1]], col = cols[1], lty = ltys[1])
for (i in 2:length(Coverages)) {
    coverage.uniformity(Coverages[[i]], add = T, col = cols[i], lty = ltys[i])
}
legend('bottomleft', legend = names(Coverages), lty = ltys, col = cols)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)
cat(" done\n")

# bait GC Content #
cat('Calculating bait GC content ')
baitGcContent <- NULL
for (i in 1:nrow(baits)) {
    Start <- start(baits)[i]
    End <- end(baits)[i]
    chr <- as.character(baits[["space"]])[i]
    chrName <- paste('chr', sub('chr', '', chr), sep = '')

    chrString <- DNAString(Hsapiens[[chrName]])

    baitGcContent <- c(baitGcContent, sum(letterFrequency(chrString[Start:End], c("G", "C"))) / (End - Start))
}
cat(' done\n')

# targets GC Content #
cat('Calculating targets GC content ')
targetGcContent <- NULL
for (i in 1:nrow(targets)) {
    Start <- start(targets)[i]
    End <- end(targets)[i]
    chr <- as.character(targets[["space"]])[i]
    chrName <- paste('chr', sub('chr', '', chr), sep = '')

    chrString <- DNAString(Hsapiens[[chrName]])

    targetGcContent <- c(targetGcContent, sum(letterFrequency(chrString[Start:End], c("G", "C"))) / (End - Start))
}
cat(' done\n')


cat("Generating coverage strip plots\n")

cat("\tacross baits...\n")
bcn <- do.call(cbind, baitcov.norms)
rownames(bcn) <- paste(baits[["space"]], ":", start(baits), "-", end(baits), sep = '')
oo <- order(as.integer(sub('chr', '', as.character(baits[["space"]]))), as.integer(start(baits)), decreasing = F)
bcn <- bcn[oo, ]
baitGcContent <- baitGcContent[oo]
bclist <- lapply(apply(bcn, 1, list), unlist)
gfn <- paste(opt$outDir, '/baitcov.png', sep = '')
png(gfn, height = 1000, width = 30 * nrow(bcn), type = 'cairo-png')
par(mar = c(15,5,5,5))
beeswarm(bclist, horizontal = F, las = 2, ylab = "Normalized Coverage", pch = 1, ylim = c(0, quantile(bcn,0.95)))
par(new = T)
beeswarm(as.list(baitGcContent), yaxt = 'n', xaxt = 'n', col = 'red', pch = 3, ylim = c(0, 1))
axis(4)
mtext("GC content", 4, line = 3)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

cat("\tacross targets ...\n")
tcn <- do.call(cbind, targetcov.norms)
rownames(tcn) <- paste(targets[["space"]], ":", start(targets), "-", end(targets), sep = '')
oo <- order(as.integer(sub('chr', '', as.character(targets[["space"]]))), as.integer(start(targets)), decreasing = F)
tcn <- tcn[oo, ]
targetGcContent <- targetGcContent[oo]
tclist <- lapply(apply(tcn, 1, list), unlist)
gfn <- paste(opt$outDir, '/targetcov.png', sep = '')
png(gfn, height = 1000, width = 30 * nrow(tcn), type = 'cairo-png')
par(mar = c(15,5,5,5))
beeswarm(tclist, horizontal = F, las = 2, ylab = 'Normalized Coverage', ylim = c(0, quantile(tcn, 0.95)))
par(new = T)
beeswarm(as.list(baitGcContent), yaxt = 'n', xaxt = 'n', col = 'red', pch = 3, ylim = c(0, 1))
axis(4)
mtext("GC content", 4, line = 3)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)


#plot(bcn[, 1], 1:nrow(bcn), yaxt = 'n', type = 'b', ylab = '', xlab = 'Normalized Coverage')
#for (i in 2:ncol(bcn)) {
    #points(bcn[, i], 1:nrow(bcn), col = i, type = 'b')
#}
#lines(baitGcContent, 1:nrow(bcn), lty = 3)
#axis(2, at = 1:nrow(bcn), labels = rownames(bcn), las = 2)

cat("\tper sample ...\n")
bclist2 <- lapply(apply(bcn, 2, list), unlist)
gfn <- paste(opt$outDir, '/samplebaitcov.png', sep = '')
png(gfn, height = 100 * ncol(bcn), width = 1000, type = 'cairo-png')
par(mar = c(5,15,5,5))
boxplot(bclist2, horizontal = T, outline = F, las = 2, col = 'red')
beeswarm(bclist2, horizontal = T, add = T)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

cat('done\n')

cat("Plotting coverage vs GC ... ")

gfn <- paste(opt$outDir, '/baitGCcontent.png', sep = '')
png(gfn, height = 800, width = 800, type = 'cairo-png')
plot(baitGcContent, bcn[,1], xlab = "GC Content", ylab = "Normalized Coverage", pch = 19, ylim = c(0, max(bcn)))
for (i in 2:ncol(bcn)) {
    points(baitGcContent, bcn[,i], pch = 19)
}
abline(h = 1, lty = 2)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)


cat(' done\n')

cat("Plotting GC correlation ")

gcCor <- apply(bcn, 2, function(x) cor(x, baitGcContent, method = 'spear'))
gfn <- paste(opt$outDir, '/gcCorBarplot.png', sep = '')
png(gfn, height = 1200, width = 1200, type = 'cairo-png')
par(oma = c(5,15,5,5))
barplot(gcCor, horiz = T, las = 1, xlab = "Spearman Correlation Coefficient", main = "GC Content Correlation")
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

cat(" done\n")

cat("Plotting coverage correlation heatmap ... ")

gfn <- paste(opt$outDir, '/covCorHeatmap.png', sep = '')
png(gfn, height = 1200, width = 1200, type = 'cairo-png')
cols <- brewer.pal(9, 'Blues')
heatmap.2(cor(bcn), trace = 'none', scale = 'none', mar = c(10, 10), col = cols, main = "Sample Coverage Correlation")
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

cat(" done\n")

cat("Generating coverage histograms and GC correlation plots ")
for (n in names(Coverages)) {
    gfn <- paste(opt$outDir, '/', n, '.covHist.png', sep = '')
    png(gfn, height = 800, width = 800, type = 'cairo-png')
    if (!is.null(opt$covThreshold)) {
        coverage.hist(Coverages[[n]]$coverageTarget, covthreshold = opt$covThreshold, main = n)
    } else {
        coverage.hist(Coverages[[n]]$coverageTarget, main = n)
    }
    null <- dev.off()
    hwriteImage(basename(gfn), pg)

    gfn <- paste(opt$outDir, '/', n, '.baitGCcontent.png', sep = '')
    png(gfn, height = 800, width = 800, type = 'cairo-png')
    plot(baitGcContent, bcn[,n], xlab = "GC Content", ylab = "Normalized Coverage", pch = 19, ylim = c(0, max(bcn)), main = n)
    lines(lowess(baitGcContent, bcn[, n]))
    null <- dev.off()
    hwriteImage(basename(gfn), pg, br = T)

    cat('.')
}

cat(' done\n')

closePage(pg)
