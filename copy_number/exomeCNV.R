#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("doMPI"))
suppressPackageStartupMessages(library("multicore"))
suppressPackageStartupMessages(library("ExomeCNV"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, type = "character", action = "store", help ="Output directory (required)"),
                make_option("--useMPI", default = F, action = "store_true", help ="Use MPI [default %default]"),
                make_option("--numThreads", default = 1, type = "integer", action = "store", help ="Number of threads [default %default]"),
                make_option("--admixtureRate", default = 0.5, type = "double", action = "store", help ="Sample admixture rate [default %default]"),
                make_option("--sensSpec", default = 0.9999, type = "double", action = "store", help ="Minimum sensitivity and specificity to call an exon CNV [default %default]"),
                make_option("--cbsSensSpec", default = 0.99, type = "double", action = "store", help ="CBS sensitivity and specificity [default %default]"),
                make_option("--readLen", default = 75, type = "integer", action = "store", help ="Number of threads [default %default]"))

parser <- OptionParser(usage = "%prog [options] tumorGatkCoverageFile normalGatkCoverageFile", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need GATK coverage file (sample_interval_summary) for tumor and normal\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n\n");
    print_help(parser);
    stop();
} else {
    tumorFile <- arguments$args[1];
    normalFile <- arguments$args[2];
}

if (opt$useMPI) {
    cl <- startMPIcluster()
    registeDoMPI(cl)
} else {
    registerDoMC(opt$numThreads)
}

read.coverage.gatk <- function (file) {
    gatk = read.table(file, header = TRUE)
    gatk <- gatk[grepl('-', as.character(gatk$Target)), ]
    chrpos = matrix(unlist(strsplit(as.character(gatk$Target),
                                    ":")), ncol = 2, byrow = TRUE)
    chr = factor(paste("chr", chrpos[, 1], sep = ""))
    pos = matrix(as.integer(unlist(strsplit(chrpos[, 2], "-"))),
                 ncol = 2, byrow = TRUE)
    start = pos[, 1]
    end = pos[, 2]
    return(data.frame(probe = gatk$Target, chr = chr, probe_start = start,
                      probe_end = end, targeted.base = end - start + 1, sequenced.base = NA,
                      coverage = as.numeric(gatk$total_coverage), average.coverage = as.numeric(gatk$average_coverage),
                      base.with..10.coverage = NA))
}


normalName <- sub('\\..*', '', sub('.*/', '', normalFile))
tumorName <- sub('\\..*', '', sub('.*/', '', tumorFile))
normal <- read.coverage.gatk(normalFile)
tumor <- read.coverage.gatk(tumorFile)

chrs <- levels(normal$chr)

# analysis overview
#1. Calculate log coverage ratio between case and control
#2. Call CNV for each exon individually
#3. Combine exonic CNV into segments using Circular Binary Segmentation (CBS)
#4. plot/export results

# log coverage
logR <- calculate.logR(normal, tumor)

# Call CNV for each exon individually
cnv <- foreach(i = 1:length(chrs), .combine = rbind) %dopar% {
    idx <- (normal$chr == chrs[i]);
    classify.eCNV(normal = normal[idx,], tumor = tumor[idx,], logR = logR[idx], min.spec = opt$sensSpec, min.sens = opt$sensSpec, option="spec", c = opt$admixtureRate, l = opt$readLen);
}

# Combine exonic CNV into segments using Circular Binary Segmentation (CBS)
eCNV = multi.CNV.analyze(normal, tumor, logR=logR, all.cnv.ls=list(cnv), coverage.cutoff=5, min.spec=opt$cbsSensSpec, min.sens=opt$cbsSensSpec, option="auc", c = opt$admixtureRate)

# plot and export results
name <- paste(opt$outDir, "/", tumorName, "_", normalName, sep = "")
`%+%` <- function(x, y) paste(x, y, sep = "")
write.table(cnv, file = name %+% ".cnv.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)
write.table(eCNV[!is.nan(eCNV$logR) & !eCNV$logR %in% c(-Inf,
                                                        Inf), c("chr", "probe_start", "probe_end", "logR")],
            file = name %+% ".exon.lrr.txt", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(cnv[!is.nan(cnv$logR) & !cnv$logR %in% c(-Inf,
                                                     Inf), c("chr", "probe_start", "probe_end", "logR")],
            file = name %+% ".segment.lrr.txt", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
write.table(cnv[!is.nan(cnv$logR) & !cnv$logR %in% c(-Inf, Inf) & !is.na(cnv$copy.number), c("chr", "probe_start", "probe_end", "copy.number")], file = name %+% ".segment.copynumber.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
png(filename = name %+% ".cnv.png", res = 70, width = 2000,
    height = 1200, pointsize = 16, type = "cairo-png")
do.plot.eCNV(cnv, style = "bp", lim.quantile = 0.999, bg.cnv = data.frame(chr = eCNV$chr,
                                                                          probe_end = (eCNV$probe_end + eCNV$probe_start)/2, logR = eCNV$logR),
             line.plot = TRUE)
dev.off()

if (opt$useMPI) {
    closeCluster(cl)
}
