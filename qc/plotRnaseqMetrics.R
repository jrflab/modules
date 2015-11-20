#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("hwriter"));
suppressPackageStartupMessages(library("RColorBrewer"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }));
}

optList <- list(
                make_option("--outDir", default = ".", help = "Output dir"));

parser <- OptionParser(usage = "%prog [options] [rnaseq_metrics] [normalized_coverage.rnaseq_metrics]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;


if (length(arguments$args) != 2) {
    cat("Need input rnaseq metric files\n");
    print_help(parser);
    stop();
} else {
    rsMetricsFile <- arguments$args[1];
    normCovMetricsFile <- arguments$args[2];
}

normCovMetrics <- read.table(normCovMetricsFile, sep = '\t', header = T, as.is = T)

pg <- openPage("index.html", dirname = opt$outDir)

gfn <- paste(opt$outDir, "/normalizedCoverageHistograms.png", sep = "")
png(gfn, height = 1200, width = 1200, type = 'cairo-png')
n <- ceiling((ncol(normCovMetrics) - 1) / 8)
ltys <- as.integer(gl(n, 8))
cols <- rep(1:8, n)
plot(normCovMetrics[, 1], normCovMetrics[, 2], type = 'l', col = cols[1], lty = ltys[1])
if (ncol(normCovMetrics) > 2) {
    for (i in 3:(ncol(normCovMetrics)-1)) {
        lines(normCovMetrics[, 1], normCovMetrics[, i], type = 'l', col = cols[i-1], lty = ltys[i-1])
    }
}
legend('bottom', legend = sub('\\..*', '', colnames(normCovMetrics)[-1], perl = T), lty = ltys, col = cols)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)


rsMetrics <- read.table(rsMetricsFile, sep = '\t', header = T, fill = T, row.names = 'SAMPLE', as.is = T)
colsToPlot <- c("PF_BASES", "PF_ALIGNED_BASES", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "PCT_USABLE_BASES", "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS")


gfn <- paste(opt$outDir, "/region_pcnt_barplot.png", sep = "")
png(gfn, height = 300 + 20 * nrow(rsMetrics), width = 600, type = 'cairo-png')
par(mar = c(5,10,5,5), oma = c(0,0,0,0))
stack <- c("PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES")
legtext <- c("% Ribosomal", "% Coding", "% UTR", "% Intronic", "% Intergenic")
oo <- order(rsMetrics[,"PCT_CODING_BASES"])
cols <- brewer.pal(5, 'Set1')
barplot(t(rsMetrics[oo, stack]), horiz = T, las = 2, legend.text = legtext, col = cols, args.legend = list(x = 'top', horiz = T, inset = c(-0.1, -0.1)))
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

for (cp in colsToPlot) {
    gfn <- paste(opt$outDir, "/", tolower(cp), "_barplot.png", sep = "")
    png(gfn, height = 300 + 20 * nrow(rsMetrics), width = 600, type = 'cairo-png')
    par(mar = c(5,10,5,5))
    barplot(rsMetrics[,cp], names.arg = rownames(rsMetrics), horiz = T, las = 2, main = cp)
    null <- dev.off()
    hwriteImage(basename(gfn), pg, br = T)
}
closePage(pg)

