#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("hwriter"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = ".", help = "Output dir"))

parser <- OptionParser(usage = "%prog [options] [hs_metrics.txt]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input hs metrics file\n");
    print_help(parser);
    stop();
} else {
    f <- arguments$args[1];
}

hsMetrics <- read.table(f, header = T, row.names = 1, sep = '\t')
colsToPlot <- c("TOTAL_READS", "PCT_PF_UQ_READS", "PCT_PF_UQ_READS_ALIGNED", "PCT_OFF_BAIT", "MEAN_TARGET_COVERAGE", "FOLD_ENRICHMENT", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_50X", "AT_DROPOUT", "GC_DROPOUT")

for (cp in colsToPlot) {
    gfn <- paste(opt$outDir, "/", tolower(cp), "_barplot.pdf", sep = "")
    pdf(gfn, height = 3 + nrow(hsMetrics) / 2, width = 6)
    par(mar = c(5,10,5,5))
    barplot(hsMetrics[,cp], names.arg = rownames(hsMetrics), horiz = T, las = 2, main = cp)
    null <- dev.off()
}

