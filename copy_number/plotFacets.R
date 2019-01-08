#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("facets"))
suppressPackageStartupMessages(library("pctGCdata"))

source('modules/copy_number/facetsPlotSampleLRR.R')

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}


optList <- list(
                make_option("--centromereFile", default = NULL, type = "character", action = "store", help ="centromere file"),
                make_option("--pqLine", default = F, action = "store_true", help = "draw pq centromere line"),
                make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [facets Rdata file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need facets Rdata file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    facetsFile <- arguments$args[1];
}

load(facetsFile)

tumorName <- facetsFile %>% sub('.*/', '', .) %>% sub('_.*', '', .)
normalName <- facetsFile %>% sub('.*/', '', .) %>% sub('^.*_', '', .) %>% sub('\\..*', '', .)

pdf(file = str_c(opt$outPrefix, ".lrr.pdf"), height = 3, width = 9)
plotSampleLRR(out2, fit)
dev.off()

pdf(file = str_c(opt$outPrefix, ".cncf.pdf"), height = 9, width = 9)
plotSample(out2, fit)
dev.off()

df <- left_join(out2$jointseg, out2$out)
df$chrom <- as.character(df$chrom)
df$chrom[df$chrom == "23"]  <- "X"
df$chrom[df$chrom == "24"]  <- "Y"

colours <- df$tcn
colours[df$tcn == 2] <- "black"
colours[df$tcn == 1] <- "darkred"
colours[df$tcn == 3] <- "darkgreen"
colours[df$tcn >= 4] <- "green"
colours[df$tcn == 0] <- "red"

ylim <- c(min(df$cnlr), max(df$cnlr) + 0.5)
ylim[2] <- ylim[2]+0.5
pdf(paste(opt$outPrefix,".cnlr_plot.pdf", sep=""), height=5, width=18)
plot(df$cnlr, pch=20, xlab='Index', ylab = "Copy number log-ratio", col = colours, ylim = ylim)
abline(v=cumsum(rle(df$chrom)$lengths), col = "red", lty = 3)
text(cumsum(rle(df$chrom)$lengths)-((rle(df$chrom)$lengths)/2), ylim[2]-0.25, labels = unique(df$chrom))

if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1] == j)[1],3]
        index <- which(df$chrom == j & df$maploc > pos)[1]
        if (opt$pqLine && !is.na(index)) {
            abline(v = index, col = "darkgrey", lty = 3)
        }
    }
}
dev.off()


for (chr in unique(df$chrom)) {
    chdf <- filter(df, chrom == chr)
    colours <- chdf$tcn
    colours[chdf$tcn == 2] <- "black"
    colours[chdf$tcn == 1] <- "darkred"
    colours[chdf$tcn == 3] <- "darkgreen"
    colours[chdf$tcn >= 4] <- "green"
    colours[chdf$tcn == 0] <- "red"

    if (nrow(chdf) > 0) {
        ylim <- c(min(chdf$cnlr), max(chdf$cnlr) + 0.5)
        pdf(paste(opt$outPrefix,".cnlr_plot.", chr, '.pdf', sep=""), height = 5, width = 6)
        plot(chdf$cnlr, pch=20, xlab='Index', ylab="Copy number", ylim=ylim, col = colours, main = paste('Chromosome', chr))
        points(chdf$cnlr.median.clust, pch = 20, col = 'blue')

        if (!is.null(opt$centromereFile)) {
            cen <- read.table(opt$centromereFile, sep = '\t')
            for (j in unique(cen[,1])) {
                pos <- cen[which(cen[,1]==j)[1],3]
                index <- which(chdf$chrom == j & chdf$maploc > pos)[1]
                if (opt$pqLine && !is.na(index)) {
                    abline(v = index, col = "darkgrey", lty = 3)
                }
            }
        }
        dev.off()
    }
}

