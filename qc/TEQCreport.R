#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("TEQC"))
suppressPackageStartupMessages(library("hwriter"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomeGraphs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--offset", default = 100, help ="target offset [default %default]"),
                make_option("--plotOffset", default = 1000, help ="coverage plot offset [default %default]"),
                make_option("--covThreshold", default = NULL, help = "coverage hist threshold [default %default]"),
                make_option("--variantPosFile", default = NULL, help = "tab-delimited file containing variant positions (chr pos)[default %default]"),
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

if (!is.null(opt$variantPosFile)) {
    cat("Loading variant positions ... ")
    variantPos <- read.table(opt$variantPosFile, sep = '\t', header = F, as.is = T, col.names = c("chr", "pos"))
    cat(" done\n")
}


avgcovs <- list()
Coverages <- list()
enrichments <- list()
frs <- list()
targetcovs <- list()

cat("Loading R data: ")
for (f in files) {
    n <- sub('\\..*', '', basename(f))
    cat(n, " ")
    load(f)
    avgcovs[[n]] <- avgcov
    Coverages[[n]] <- Coverage
    enrichments[[n]] <- enrichment
    targetcovs[[n]] <- targetcov
    frs[[n]] <- fr
}
cat(' done\n')

avgcov <- do.call('rbind', avgcovs)

pg <- openPage('index.html', dirname = opt$outDir, title = 'QC report')

write.table(avgcov, file = paste(opt$outDir, '/avgcov.txt', sep = ''), quote = F, sep = '\t')
hwrite(avgcov, page = pg, br = T)

cat("Generating coverage uniformity plots ... ")
gfn <- paste(opt$outDir, '/coverageUniformity.png', sep = '')
png(gfn, height = 600, width = 600)
coverage.uniformity(Coverages[[1]])
for (i in 2:length(Coverages)) {
    coverage.uniformity(Coverages[[i]], add = T, col = i)
}
legend('bottomleft', legend = names(Coverages), lty = 1, col = 1:length(Coverages))
dev.off()
hwriteImage(basename(gfn), pg, br = T)
cat(" done\n")

cat("Generating coverage histograms ")
for (n in names(Coverages)) {
    gfn <- paste(opt$outDir, '/', n, '.covHist.png', sep = '')
    png(gfn, height = 600, width = 600)
    if (!is.null(opt$covThreshold)) {
        coverage.hist(Coverages[[n]]$coverageTarget, covthreshold = opt$covThreshold, main = n)
    } else {
        coverage.hist(Coverages[[n]]$coverageTarget, main = n)
    }
    dev.off()
    hwriteImage(basename(gfn), pg, br = T)
    cat('.')
}
cat(' done\n')

cat("Generating coverage plots per bait ")
windowSize <- 1000
for (i in 1:nrow(baits)) {
    Start <- start(baits)[i]
    End <- end(baits)[i]
    chr <- as.character(space(baits)[i])
    chrName <- paste('chr', sub('chr', '', chr), sep = '')
    chrString <- DNAString(Hsapiens[[chrName]])
    gcContent <- rowSums(letterFrequencyInSlidingView(chrString[Start:(End + windowSize - 1)], windowSize, c("G", "C"))) / windowSize

    ir <- IRanges(start = Start, end = End)
    xlab <- paste("Chromosome", chr)
    ylab <- "Normalized Coverage"

    gfn <- paste(opt$outDir, '/coveragePlot.', i, '.png', sep = '')
    png(gfn, height = 400 * length(Coverages), width = 2000)
    par(mfrow = c(length(Coverages), 1), mar = c(0, 5, 0, 5), oma = c(5,5,5,5))
    for (j in 1:length(Coverages)) {
        covercounts <- Coverages[[j]]$coverageAll[[chr]]
        covsel <- seqselect(covercounts, ir)
        covsel <- covsel / mean(covsel)
        plot(Start:End, covsel, type = "l",
             xlab = xlab, ylab = ylab, axes = F, ylim = c(0,3), col = 1, lty = 1)
        #lines(lowess(Start:End, covsel, f = 0.05), col = 'blue', lty = 3,  lwd = 3)
        axis(side = 2, at = pretty(range(0:3)))
        mtext(names(Coverages)[j], side = 2, line = 5)
        abline(h = 1, lty = 3)
        if (!is.null(opt$variantPosFile)) {
            # plot variant positions
            posns <- subset(variantPos, chr == chr & pos > Start & pos < End)
            for (pos in posns) {
                abline(v = pos, col = 'blue', lty = 2)
            }
        }
        par(new = T)
        plot(Start:End, gcContent, axes = F, type = 'l', xlab = '', ylab = '', bty = 'n', col = 2, lty = 2)
        axis(side = 4, at = pretty(range(gcContent)))
        mtext("GC", side = 4, line = 3)
    }
    axis(side = 1, at = pretty(range(Start:End)))
    mtext(xlab, side = 1, line = 3)
    dev.off()
    hwriteImage(basename(gfn), pg)

    gfn <- paste(opt$outDir, '/coverageGCPlot.', i, '.png', sep = '')
    png(gfn, height = 400 * length(Coverages), width = 600)
    par(mfrow = c(length(Coverages), 1), mar = c(0,5,0,5), oma = c(5,5,5,5))
    for (j in 1:length(Coverages)) {
        covercounts <- Coverages[[j]]$coverageAll[[chr]]
        covsel <- seqselect(covercounts, ir)
        smoothScatter(gcContent, covsel, ylab = "Coverage", axes = F)
        lines(lowess(gcContent, covsel), lty = 2)
        axis(side = 2, at = pretty(range(covsel)))
        mtext(names(Coverages)[j], side = 2, line = 5)
    }
    axis(side = 1, at = pretty(range(gcContent)))
    mtext("GC", side = 1, line = 3)
    dev.off()
    hwriteImage(basename(gfn), pg, br = T)

    cat(".")
}
cat(' done\n')

closePage(pg)
