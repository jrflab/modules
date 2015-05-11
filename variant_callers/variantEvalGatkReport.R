#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("gsalib"));
suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("lattice"));
suppressPackageStartupMessages(library("hwriter"));
suppressPackageStartupMessages(library("maptools"));
suppressPackageStartupMessages(library("corrgram"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }));

optList <- list(
                make_option("--width", default = 1000, action = 'store', help = "Plot width [default %default]"),
                make_option("--height", default = 1000, action = 'store', help = "Plot height [default %default]"),
                make_option("--metrics", default = NULL, type = 'character', action = 'store', help = "Picard metrics matrix"),
                make_option("--outDir", default = NULL, type = "character", action = "store", help ="Output directory (required)"))

parser <- OptionParser(usage = "%prog [options] variantEvalFile", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need GATK variant eval file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$metrics)) {
    cat("Need metrics matrix (bamIntervalMetrics.mk:hs_metrics.txt)\n\n");
    print_help(parser);
    stop();
} else {
    variantEvalFile <- arguments$args[1];
}

pg <- openPage('index.html', dirname = opt$outDir)

report <- gsa.read.gatkreport(variantEvalFile);

metricsData <- read.table(opt$metrics, header = T, row.names = 1, as.is = T, sep = '\t')



# comp overlap
# dbsnp known
#compOverlap <- subset(report[["CompOverlap"]], Novelty != "all" & IntervalStratification == "overlaps.intervals" & Filter == 'called')
#compOverlap$Novelty <- factor(compOverlap$Novelty)
#compOverlap$IntervalStratification <- factor(compOverlap$IntervalStratification)
compOverlap <- subset(report[["CompOverlap"]], Novelty != "all" & Filter == 'called')
compOverlap$Novelty <- factor(compOverlap$Novelty)

#barchart(nEvalVariants ~ EvalRod | IntervalStratification, group = Novelty, data = compOverlap, scales = list(x = list(rot = 90)), origin = 0, auto.key = list(space = 'top', rectangles = T))

gfn <- paste(opt$outDir, "/knownConcordancy.png", sep = "")
trellis.device(device = "png", filename = gfn, height = opt$height, width = opt$width)
#trellis.par.set(fontsize = list(text = 12))
lb <- min(subset(compOverlap, Novelty == 'known')$concordantRate)
bc <- barchart(EvalRod ~ concordantRate, data = subset(compOverlap, Novelty == 'known'), origin = 0, xlim = c(lb-5, 100),  main = "Concordancy at Known Sites", xlab = "Concordant Rate")
print(bc)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

#countVariants <- subset(report[["CountVariants"]], Novelty != 'all' & IntervalStratification == 'overlaps.intervals' & Filter == 'called')
#countVariants$IntervalStratification <- factor(countVariants$IntervalStratification)
countVariants <- subset(report[["CountVariants"]], Novelty != 'all' & Filter == 'called')
countVariants$Novelty <- factor(countVariants$Novelty)



gfn <- paste(opt$outDir, "/knownConcordancyVsNumVar.png", sep = "")
trellis.device(device = "png", filename = gfn, height = opt$height, width = opt$width)
snpConcordantRate <- subset(compOverlap, Novelty == 'known')$concordantRate
names(snpConcordantRate) <- subset(compOverlap, Novelty == 'known')$EvalRod
nVariants <- subset(countVariants, Novelty == 'novel')$nVariantLoci
names(nVariants) <- subset(countVariants, Novelty == 'novel')$EvalRod
xl <- range(snpConcordantRate) + c(-5, 5)
labs <- as.character(subset(compOverlap, Novelty == 'known')$EvalRod)
p <- xyplot(nVariants ~ snpConcordantRate, xlim = xl, pch = 19, labels = labs, panel = function(x, y, labels, ...) {
       panel.xyplot(x, y, ...)
       panel.pointLabel(x, y, labels = labels, ...)
})
print(p)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

gfn <- paste(opt$outDir, "/metricsCorrgram.png", sep = "")
metrics <- c("SNP_CONCORDANT_RATE", "NUM_VARIANTS", "GC_DROPOUT", "PCT_OFF_BAIT", "MEAN_TARGET_COVERAGE", "PCT_USABLE_BASES_ON_BAIT")
trellis.device(device = "png", filename = gfn, height = 400 * length(metrics), width = 400 * length(metrics))
metricsData <- cbind(metricsData, SNP_CONCORDANT_RATE = snpConcordantRate[rownames(metricsData)], NUM_VARIANTS = nVariants[rownames(metricsData)]) 
corrgram(metricsData[,metrics], order = T, lower.panel = panel.ellipse, upper.panel = panel.pts, text.panel = panel.txt, diag.panel = panel.minmax)
        hwriteImage(basename(gfn), pg, br = T)
null <- dev.off()


for (metric in colnames(countVariants)[9:ncol(countVariants)]) {
    if (sum(countVariants[, metric]) > 0) {
        gfn <- paste(opt$outDir, "/", metric, '.barchart.png', sep = '')
        trellis.device(device = "png", filename = gfn, height = opt$height, width = opt$width)
        form <- as.formula(paste('EvalRod ~', metric))
        bc <- barchart(form, groups = Novelty, auto.key = list(space = "top", rectangles = T), data = countVariants, origin = 0)
        print(bc)
        null <- dev.off()
        hwriteImage(basename(gfn), pg, br = T)
    }
}

## TiTv ##
TiTv <- subset(report[["TiTvVariantEvaluator"]], Novelty != 'all' & Filter == 'called')
TiTv$Novelty <- factor(TiTv$Novelty)
metrics <- c("nTi", "nTv", "tiTvRatio")
for (metric in metrics) {
    gfn <- paste(opt$outDir, "/", metric, '.barchart.png', sep = '')
    trellis.device(device = "png", filename = gfn, height = opt$height, width = opt$width)
    form <- as.formula(paste('EvalRod ~', metric))
    bc <- barchart(form, groups = Novelty, auto.key = list(space = "top", rectangles = T), data = TiTv, origin = 0)
    print(bc)
    null <- dev.off()
    hwriteImage(basename(gfn), pg, br = T)
}

## Indel length hists ##
gfn <- paste(opt$outDir, "/indelLengthHists.png", sep = '')
indelLength <- subset(report[["IndelLengthHistogram"]], Novelty == 'all' & Filter == 'called')
n <- ceiling(sqrt(nlevels(indelLength$EvalRod)))
trellis.device(device = "png", filename = gfn, height = opt$height * n, width = opt$width * n)
xyp <- xyplot(Freq ~ Length | EvalRod, data = indelLength, type = 'b')
print(xyp)
null <- dev.off()
hwriteImage(basename(gfn), pg, br = T)

## Indel summary ##
indelSummary <- subset(report[["IndelSummary"]], Novelty != 'all' & Filter == 'called')
indelSummary$Novelty <- factor(indelSummary$Novelty)
metrics <- c("n_indels", "n_singleton_indels", "SNP_to_indel_ratio", "SNP_to_indel_ratio_for_singletons", "n_insertions", "n_deletions", "insertion_to_deletion_ratio") 
for (metric in metrics) {
    gfn <- paste(opt$outDir, "/", metric, '.barchart.png', sep = '')
    trellis.device(device = "png", filename = gfn, height = opt$height, width = opt$width)
    form <- as.formula(paste('EvalRod ~', metric))
    bc <- barchart(form, groups = Novelty, auto.key = list(space = "top", rectangles = T), data = indelSummary, origin = 0)
    print(bc)
    null <- dev.off()
    hwriteImage(basename(gfn), pg, br = T)
}


closePage(pg)
