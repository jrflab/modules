#!/usr/bin/env Rscript
# plots control freec 

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--out", default = NULL, help = "Output file pdf"),
                make_option("--groupByGene", action = 'store_true', default = F, help = "do not collapse by cytoband"),
                make_option("--centromereMatrix", help = "Centromere position matrix"));

parser <- OptionParser(usage = "%prog [options] allthresolded.by_genes.txt", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input gistic allthresholded.by_genes.txt file\n");
    print_help(parser);
    stop();
}
if (is.null(opt$out)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
}



collapseByCytoband <- function(tab, include_subbands=T) {
    if (!include_subbands) {
        tab$Cytoband <- unlist(lapply(tab$Cytoband, function(x) { strsplit(x, split=".", fixed=T)[[1]][1]}))}
    cyto <- unique(tab$Cytoband)    

    res <- matrix(nrow=length(unique(tab$Cytoband)), ncol=0)
    for (i in 3:ncol(tab)) {
        res <- cbind(res, tapply(tab[,i], tab$Cytoband, function(x) { round(median(x))}))
        colnames(res)[ncol(res)] <- colnames(tab)[i]
    }

    res <- res[match(cyto, rownames(res)),]
    res
}

plotHeatmap <- function(all_thresholded_file, plotfile, pheno=NULL, genes=NULL, cytoband=NULL, group.by="cytoband", ...) {

    mat <- read.delim(all_thresholded_file, as.is=T, row.names=1)
    if (!is.null(genes)) {mat <- mat[which(row.names(mat) %in% genes),]}


    if (group.by=="cytoband") {
        mat2 <- collapseByCytoband(mat, ...)
        if (!is.null(cytoband)) {
            mat2 <- mat2[which(rownames(mat2) %in% cytoband),]}
        chr <- unlist(lapply(rownames(mat2), function(x) {strsplit(x, split="p|q", perl=T)[[1]][1]}))
        } else { mat2 <- as.matrix(mat[,-c(1,2)])
        chr <- unlist(lapply(mat$Cytoband, function(x) {strsplit(x, split="p|q", perl=T)[[1]][1]}))
    }

    chrsep <- cumsum(rle(chr)$lengths)
    chrmid <- c(0,chrsep[-length(chrsep)]) + (rle(chr)$lengths/2)


    if (!is.null(pheno)) {mat2 <- mat2[,match(pheno, colnames(mat2))]
    } else { mat2 <- mat2[,order(colnames(mat2))]}
    mat2 <- mat2[,ncol(mat2):1]

    pdf(plotfile, height=max(3, ncol(mat2)/3), width=12)
    par(mar=c(6,10,1,2))
    image(mat2, col=c("red", "darksalmon", "white", "lightblue", "blue"), xaxt='n', yaxt='n', zlim=c(-2,2))
    box()

    if (!is.null(cytoband)) {
        axis(1,at=seq(0,1,1/(nrow(mat2)-1)), label=rownames(mat2), las=2, cex.axis=1, tick=F)
    } else {
        axis(1,at=chrmid/(max(chrsep)-1), label=rle(chr)$values, cex.axis=0.8, tick=F)
    }
    axis(2,at=seq(0,1,1/(ncol(mat2)-1)), label=colnames(mat2), las=2, cex.axis=1, tick=F)

    for (i in (chrsep*2)-1) {
        abline(v=i/((max(chrsep)-1)*2), col="grey")
    }

    for (i in seq(-1, ((2*(ncol(mat2)-1))+1), 2)) {
        abline(h=i/(2*(ncol(mat2)-1)), col="white", lwd=2)
    }

    null <- dev.off()
}



groupBy <- ifelse(opt$groupByGene, "gene", "cytoband")
fn <- arguments$args[1];
plotHeatmap(fn, opt$out, group.by = groupBy)

