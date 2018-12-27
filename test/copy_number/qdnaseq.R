#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("QDNAseq"))
suppressPackageStartupMessages(library("future"))

future::plan("multiprocess")
options(mc.cores=16L)


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

if (!dir.exists("qdnaseq/readcounts")) {
	dir.create("qdnaseq/readcounts")
}

if (!dir.exists("qdnaseq/isobars")) {
	dir.create("qdnaseq/isobars")
}

if (!dir.exists("qdnaseq/variance")) {
	dir.create("qdnaseq/variance")
}

if (!dir.exists("qdnaseq/log2ratio")) {
	dir.create("qdnaseq/log2ratio")
}

if (!dir.exists("qdnaseq/bed")) {
	dir.create("qdnaseq/bed")
}

args_list = list(make_option("--sample", default = NA, type = 'character', help = "sample name"),
				 make_option("--binsize", default = NA, type = 'character', help = "bin size"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (is.na(as.numeric(opt$binsize))) {
	opt$binsize = 30
} else {
	opt$binsize = as.numeric(opt$binsize)
}

bins = getBinAnnotations(binSize=opt$binsize, genome="hg19")
readCounts = binReadCounts(bins=bins, bamfiles=paste0("bam/", opt$sample, ".bam"),
						   isPaired=TRUE,
       					   isProperPair=TRUE,
        				   minMapq=30,
        				   pairedEnds=TRUE,
        				   chunkSize=TRUE)
       
# read counts versus genomic coordinates
pdf(file=paste0("qdnaseq/readcounts/", opt$sample, ".pdf"), width=14, height=9)
plot(readCounts, logTransform=TRUE, ylim=c(0, 20))
highlightFilters(readCounts, logTransform=TRUE, residual=TRUE, blacklist=TRUE)
dev.off()

readCountsFiltered = applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

# %GC content versus mappability
pdf(file=paste0("qdnaseq/isobars/", opt$sample, ".pdf"), width=7, height=7)
isobarPlot(readCountsFiltered)
dev.off()

readCountsFiltered = estimateCorrection(readCountsFiltered)

# noise (variance) versus bin coverage
pdf(file=paste0("qdnaseq/variance/", opt$sample, ".pdf"), width=7, height=7)
noisePlot(readCountsFiltered)
dev.off()

copyNumbers = correctBins(readCountsFiltered)
copyNumbersNormalized = normalizeBins(copyNumbers)
copyNumbersSmooth = smoothOutlierBins(copyNumbersNormalized)

# log2 ratio versus genomic coordinates
pdf(file=paste0("qdnaseq/log2ratio/", opt$sample, ".pdf"), width=14, height=9)
plot(copyNumbersSmooth, ylim=c(-4,4))
dev.off()

# write log2 ratio to file
exportBins(copyNumbersSmooth, file=paste0("qdnaseq/bed/", opt$sample, ".bed"), format="bed")
