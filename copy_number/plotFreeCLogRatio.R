#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outFile", default = NULL, help = "Output file prefix"),
                make_option("--ploidy", default = 2, help = "ploidy"),
                make_option("--centromereTable", help = "Centromere position table"));

parser <- OptionParser(usage = "%prog [options] [ratio file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input controlFreeC ratio.txt files\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}
args <- commandArgs()

ratio <- read.table(arguments$args[1], header = T)

ploidy <- opt$ploidy

png(opt$outFile, width = 1000, height = 500,
    units = "px", pointsize = 20, bg = "white", res = NA, type = 'cairo-png')

maxLevelToPlot <- 3
for (i in c(1:length(ratio$Ratio))) {
	if (ratio$Ratio[i]>maxLevelToPlot) {
		ratio$Ratio[i]=maxLevelToPlot;
	}
}

x <- as.character(ratio$Chromosome)
x[x == "X"] <- 23
if (any(x == "Y")) {
    x[x == "Y"] <- 24
}
if (any(x == "MT")) {
    x[x == "MT"] <- 25
}
x <- as.integer(x)
oo <- order(x, ratio$Start)
ratio <- ratio[oo, ]

xx <- rle(as.character(ratio$Chromosome))
os <- c(0, cumsum(as.numeric(ratio$Start[cumsum(xx$lengths)])))
names(os) <- xx$values
x <- ratio$Start + os[as.character(ratio$Chromosome)]
plot(x, ratio$Ratio * ploidy, ylim = c(0, maxLevelToPlot*ploidy), xlab = "position", ylab = "normalized copy number profile",pch = ".",col = 'Black', xaxt = 'n')
tt <- which(ratio$CopyNumber > ploidy)
points(x[tt], ratio$Ratio[tt]*ploidy,pch = ".",col = 'DarkGreen') # gain

tt <- which(ratio$Ratio == maxLevelToPlot & ratio$CopyNumber > ploidy)	
points(x[tt], ratio$Ratio[tt] * ploidy,pch = ".",col = 'green',cex=4)

tt <- which(ratio$CopyNumber < ploidy & ratio$CopyNumber!= -1)
points(x[tt],ratio$Ratio[tt] * ploidy,pch = ".", col = 'Darkred')

xx <- rle(as.character(ratio$Chromosome))
chrMid <- os[-length(os)] + diff(os) / 2
axis(1, at = chrMid, labels = names(os)[-length(os)], tick = F)
axis(1, at = os, tick = T, labels = F)
abline(v = os[-c(1,length(os))])

if (!is.null(opt$centromereTable)) {
    centromerePosns <- subset(read.table(opt$centromereTable, sep = '\t'), V8 == "centromere", as.is = T)
    centromerePosns <- centromerePosns[, 2:4]
    centromerePosns[, 1] <- sub('chr', '', centromerePosns[, 1])
    centromerePosns[centromerePosns[,1] == "X", 1] <- 23
    centromerePosns <- centromerePosns[-which(centromerePosns[,1] == "Y"), ]
    oo <- order(as.integer(centromerePosns[,1]))
    centromerePosns <- centromerePosns[oo, ]
    centromerePosns[centromerePosns[,1] == 23, 1] <- "X"

    cmPos <- centromerePosns[,2]
    names(cmPos) <- as.character(centromerePosns[,1])

    abline(v = cmPos + os[-length(os)], col = 'darkgrey', lty = 3)
}

null <- dev.off()
