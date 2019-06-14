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
suppressPackageStartupMessages(library("GAP"))

source('modules/copy_number/facetsPlotSampleLRR.R')

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}


optList <- list(make_option("--centromereFile", default = NULL, type = "character", action = "store", help ="centromere file"),
                make_option("--pqLine", default = F, action = "store_true", help = "draw pq centromere line"),
                make_option("--outPrefix", default = NULL, help = "output prefix"))
parser <- OptionParser(usage = "%prog [options] [facets Rdata file]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (length(arguments$args) < 1) {
    cat("Need facets Rdata file\n")
    print_help(parser)
    stop()
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop()
} else {
    facetsFile <- arguments$args[1]
}

load(facetsFile)

if (!dir.exists("facets/plots")) {
	dir.create("facets/plots")
}

if (!dir.exists("facets/plots/log2")) {
	dir.create("facets/plots/log2")
}

if (!dir.exists("facets/plots/cncf")) {
	dir.create("facets/plots/cncf")
}

if (!dir.exists("facets/plots/bychr")) {
	dir.create("facets/plots/bychr")
}


tumorName <- facetsFile %>%
			 sub('.*/', '', .) %>%
			 sub('_.*', '', .)
normalName <- facetsFile %>%
			  sub('.*/', '', .) %>%
			  sub('^.*_', '', .) %>%
			  sub('\\..*', '', .)

pdf(file = str_c(opt$outPrefix, ".pdf"), width=10, height=4.25)
plot_log2_(x=out2, y=fit, purity=fit$purity, ploidy=fit$ploidy, title = gsub("facets/plots/log2/", "", opt$outPrefix, fixed=TRUE))
dev.off()

pdf(file = str_c(gsub("log2", "cncf", opt$outPrefix, fixed=TRUE), ".pdf"), width=10, height=7)
plot_cncf_(out2, fit)
dev.off()

df <- left_join(out2$jointseg, out2$out)
df$chrom <- as.character(df$chrom)
df$chrom[df$chrom == "23"]  <- "X"
df$chrom[df$chrom == "24"]  <- "Y"

for (chr in unique(df$chrom)) {
     chdf <- filter(df, chrom == chr)
 
     if (nrow(chdf) > 0) {
        pdf(paste(gsub("log2", "bychr", x=opt$outPrefix, fixed=TRUE), "_", chr, '.pdf', sep=""), height = 5, width = 5)
        par(mar=c(5, 5, 4, 2)+.1)
        plot(chdf$cnlr, type="p", pch=20, col="grey80", xlab="", ylab="", main = "", ylim=c(-4,5), axes=FALSE, frame=TRUE)
        points(chdf$cnlr.median.clust, pch=20, col="red")
        abline(h=0, col="goldenrod3", lty=1, lwd=1)
        axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		axis(1, at = NULL, labels = NULL, cex.axis = 0.85, las = 1)
    	rect(xleft=1-1e10, xright=nrow(chdf)+10e3, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
		title(main = gsub("facets/plots/log2/", "", opt$outPrefix, fixed=TRUE), line=-1.5, cex.main=.75, font.main=1)
    	box(lwd=1.5)

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
