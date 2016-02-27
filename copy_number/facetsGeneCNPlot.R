#!/usr/bin/env Rscript
## plots facets geneCN file, comparing copy number amp/del between samples

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(optparse,RColorBrewer,GenomicRanges,plyr,dplyr,readr,stringr,tidyr,purrr,magrittr,crayon,foreach,Cairo,RMySQL,rtracklayer))

#--------------
# parse options
#--------------

parser <- OptionParser(usage = "%prog [geneCN file] [output_plot_file]");

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need one geneCN file and one output plot file\n")
    print_help(parser);
    stop();
} else {
    geneCN <- arguments$args[1]
    outFile <- arguments$args[2]
}

plot_heatmap <- function(facets_tab, plot_file, sample_names=NULL, col=c("red", "darksalmon", "white", "lightblue", "blue"), zlim=c(-2,2)) {
    mm <- facets_tab
    if (is.null(sample_names)) { sample_names <- list(sort(colnames(mm)[sapply(colnames(mm), function(x) {grepl("_EM$", x)})])) }
    print(sample_names)
    chrsep <- cumsum(rle(mm$chrom)$lengths)
    chrmid <- c(0,chrsep[-length(chrsep)]) + (rle(mm$chrom)$lengths/2)
    pdf(plot_file, width=12, height=5*length(sample_names))
    par(mfrow=c(length(sample_names),1), mar=c(8,5,1,2))
    lapply(sample_names, function(x, mm) {
        mm2 <- mm[,rev(x)]; #for (i in 1:ncol(mm2)) { mm2[,i] <- as.numeric(mm2[,i]) }
        print(colnames(mm2))
        image(as.matrix(mm2), col=col, xaxt='n', yaxt='n', zlim=zlim)
        box()
        for (i in (chrsep*2)-1) { abline(v=i/((max(chrsep)-1)*2), col="grey") }
        for (i in seq(-1, max(((2*(ncol(mm2)-1))+1),1), 2)) { abline(h=i/(2*(ncol(mm2)-1)), col="white", lwd=2)}
    axis(1,at=chrmid/(max(chrsep)-1), label=rle(mm$chrom)$values, cex.axis=0.8, tick=F)
    axis(2,at=seq(0,1,1/max((ncol(mm2)-1),1)), label=sub("T_.*", "", colnames(mm2)), las=2, cex.axis=1, tick=F)
    }, mm)
    legend("bottom", inset=c(0,-.4), legend=c("Homozygous deletion", "Loss", "Gain", "Amplification"),
           fill=c("red", "darksalmon", "lightblue", "blue"), xpd=T, ncol=2)
    dev.off()
}

plot_heatmap(read.table(geneCN, sep="\t", header=T, stringsAsFactors=F), outFile)
