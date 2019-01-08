#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("gplots"))

options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'b37', help = "genome build [default %default]"),
                make_option("--outPrefix", default = NULL, help = "output prefix [default %default]"))

parser <- OptionParser(usage = "%prog vcf.files", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

vcfFile <- arguments$args[1]


vcf <- readVcf(vcfFile, opt$genome)
gt <- geno(vcf)$GT
ad <- geno(vcf)$AD
af <- structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
X <- matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
X[is.na(af)] <- NA
X[af > 0.15 & af < 0.95] <- 1
X[af >= 0.95] <- 2
X[!gt %in% c("0/0", "0/1", "1/1")] <- NA
#plot(hclust(dist(t(X), method = 'manhattan')))

gt <- matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))

fn <- paste(opt$outPrefix, ".clust.pdf", sep = '')
pdf(fn, height = 9, width = 15)
null <- plot(hclust(dist(t(gt)), method = 'ward'))
dev.off()

fn <- paste(opt$outPrefix, ".heatmap.pdf", sep = '')
pdf(fn, height = 30, width = 30)
null <- heatmap.2(as.matrix(dist(t(gt))), scale = 'none', trace = 'none', keysize = 0.3, cexRow = 2, cexCol = 2, margins = c(20,20))
dev.off()

