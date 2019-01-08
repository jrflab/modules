#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("hwriter"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, help = "output dir (required)"),
                make_option("--mutations", default = NULL, help = "mutations file (required)"),
                make_option("--samples", default = 'emu/mutations.txt.samples' , help = "samples file [default %default]"),
                make_option("--sampleSubset", default = 'emu/sample_pairs.txt', help = "sample subset file: list of samples to plot contribution [default %default]"),
                make_option("--inPrefix", default = NULL, help = "EMu input prefix (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outDir)) {
    cat("Need output dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$inPrefix)) {
    cat("Need EMu input prefix\n");
    print_help(parser);
    stop();
} else if (is.null(opt$mutations)) {
    cat("Need EMu mutations file\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

glob <- paste(opt$inPrefix, '*_ml_spectra.txt', sep = '')
spectraFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '*_map_activities.txt', sep = '')
activityFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '*_assigned.txt', sep = '')
assignedFiles <- Sys.glob(glob)

pg <- openPage('index.html', dirname = opt$outDir, title = 'EMu results')

set.seed(002)
palette(sample(rainbow(30)))

for (fn in spectraFiles) {
    spectra <- read.table(fn, sep = ' ')
    spectra <- spectra[,-97] # remove empty col
    for (i in 1:nrow(spectra)) {
        ofn <- paste(opt$outDir, "/", basename(fn), sep = '')
        ofn <- sub('\\.txt$', paste("_", i, '.pdf', sep = ''), ofn)
        pdf(ofn, height = 8, width = 10)
        par(cex = 1.5)
        cols <- rep(c('LightBlue', 'Black', 'Red', 'Grey', 'Green', 'Magenta'), each = 16)
        barplot(t(spectra[i,]) * 100, beside = T, col = cols, border = cols, xaxt = 'n', main = paste("Signature", i), col.main = i, ylab = "% of mutations")
        labs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
        mtext(labs, side = 1, at = 1:6 * 16 - 7.5)
        null <- dev.off()
    }

    for (i in 1:nrow(spectra)) {
        ofn <- paste(opt$outDir, "/", basename(fn), sep = '')
        ofn <- sub('\\.txt$', paste("_", i, '.png', sep = ''), ofn)
        png(ofn, height = 500, width = 800, type = 'cairo-png')
        par(cex = 2)
        cols <- rep(c('LightBlue', 'Black', 'Red', 'Grey', 'Green', 'Magenta'), each = 16)
        barplot(t(spectra[i,]) * 100, beside = T, col = cols, border = cols, xaxt = 'n', main = paste("Signature", i), ylab = "% of mutations", col.main = i)
        labs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
        mtext(labs, side = 1, at = 1:6 * 16 - 7.5)
        null <- dev.off()
        hwriteImage(basename(ofn), pg, br = T)
    }
}

samples <- scan(opt$samples, what = 'character')
sampleSubset <- scan(opt$sampleSubset, what = 'character')

mutations <- read.table(opt$mutations, sep = ' ')
colnames(mutations) <- c('sample', 'chr', 'pos', 'snv')
x <- do.call('rbind', strsplit(as.character(mutations$sample), '\\.'))
colnames(x) <- c('sample', 'class')
mutations$sample <- as.factor(x[, 'sample'])
mutations$class <- as.factor(x[, 'class'])
mutations <- subset(mutations, sample %in% sampleSubset)

for (fn in assignedFiles) {
    assigned <- read.table(fn, sep = ' ')
    assigned <- as.matrix(assigned[,-ncol(assigned)])
    rownames(assigned) <- samples
    assigned <- assigned[sampleSubset, ]

    ofn <- paste(opt$outDir, "/", basename(fn), sep = '')
    ofn <- sub('\\.txt$', '.png', ofn)
    if (is.null(dim(assigned))) {
        s <- names(assigned)
        x <- data.frame(do.call('rbind', strsplit(s, '\\.')))
        colnames(x) <- c('sample', 'class')
        X <- do.call('rbind', split(assigned, x$sample))
        colnames(X) <- sub('.*\\.', '', colnames(X))
        tab <- table(factor(mutations$sample))
        tab <- tab[rownames(X)]
        oo <- order(tab)
        X <- X[oo, ]

        png(ofn, height = 1000, width = 1000, type = 'cairo-png')
        par(mar = c(5,10,5,5))
        barplot(t(X), beside = T, horiz = T, las = 2, legend =T, xlab = 'Number of Signature 1 Mutations', args.legend = list(x = 'bottomright'))
        null <- dev.off()
        hwriteImage(basename(ofn), pg, br = T)

        ofn <- sub('\\.png$', '.pdf', ofn)
        pdf(ofn, height = 12, width = 12)
        par(mar = c(5,10,5,5))
        barplot(t(X), beside = T, horiz = T, las = 2, legend =T, xlab = 'Number of Signature 1 Mutations', args.legend = list(x = 'bottomright'))
        null <- dev.off()
    } else {
        s <- rownames(assigned)
        x <- data.frame(do.call('rbind', strsplit(s, '\\.')))
        colnames(x) <- c('sample', 'class')
        oo <- order(x$sample)
        png(ofn, height = 1000, width = 1000, type = 'cairo-png')
        par(mar = c(5,10,5,5))
        barplot(t(assigned[oo, ]), horiz = T , las = 2, space = c(1,0), col = 1:5, xlab = "Number of mutations")
        null <- dev.off()
        hwriteImage(basename(ofn), pg, br = T)
    }
}


closePage(pg)
