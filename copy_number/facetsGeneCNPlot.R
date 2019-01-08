#!/usr/bin/env Rscript

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(library(optparse))

#--------------
# parse options
#--------------

optList <- list(
                make_option("--includeChrY", action="store_true", default=F, help="Include Chromosome Y (drop by default)"),
                make_option("--sampleColumnPostFix", default="_AWN_LRR_threshold$", help="Postfix of columns that represent samples"),
                make_option("--summaryConfig", default=NULL, help="Summary config yaml for ordering samples and replacing IDs"));
parser <- OptionParser(usage = "%prog [geneCN file] [output_plot_file]", option_list=optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need one geneCN file and one output plot file\n")
    print_help(parser)
    stop()
} else {
    geneCN <- arguments$args[1]
    outFile <- arguments$args[2]
}

# use showtext if installed for monospaced system font, useful for TCGA
# barcodes on y axis
if ("showtext" %in% rownames(installed.packages())) {
    suppressMessages(library(showtext))
    font.add("DejaVuSansMono", "DejaVuSansMono.ttf")
    showtext.auto()
    fontfamily <- "DejaVuSansMono"
} else {
    fontfamily <- "serif"
}

plot_heatmap <- function(facets_tab, plot_file, sample_column_postfix, fontfamily, sample_names=NULL, summaryConfig=NULL, col=c("blue", "lightblue", "white", "darksalmon", "red" ), zlim=c(-2,2)) {
    mm <- facets_tab
    if (is.null(sample_names)) { sample_names <- list(sort(colnames(mm)[sapply(colnames(mm), function(x) {grepl(paste(sample_column_postfix,"$",sep=""), x)})])) }
    # user order from summary_config.yaml if available
    if (!is.null(summaryConfig) && ("sample_order" %in% names(summaryConfig))) {
        sample_names[[1]] = sample_names[[1]][charmatch(summaryConfig$sample_order, sample_names[[1]])];
    }
    chrsep <- cumsum(rle(mm$chrom)$lengths)
    chrmid <- c(0,chrsep[-length(chrsep)]) + (rle(mm$chrom)$lengths/2)
    pdf(plot_file, width=20, height=5 + .8*length(sample_names[[1]]))
    par(mfrow=c(length(sample_names),1), mar=c(8,1*(max(sapply(sample_names,nchar))-nchar(sample_column_postfix)),1,2))
    lapply(sample_names, function(x, mm) {
        mm2 <- mm[,rev(x)]; #for (i in 1:ncol(mm2)) { mm2[,i] <- as.numeric(mm2[,i]) }
        image(as.matrix(mm2), col=col, xaxt='n', yaxt='n', zlim=zlim)
        box()
        for (i in (chrsep*2)-1) { abline(v=i/((max(chrsep)-1)*2), col="grey") }
        for (i in seq(-1, max(((2*(ncol(mm2)-1))+1),1), 2)) { abline(h=i/(2*(ncol(mm2)-1)), col="white", lwd=2)}
    axis(1,at=chrmid/(max(chrsep)-1), label=rle(mm$chrom)$values, cex.axis=0.8, tick=F, family=fontfamily)
    axis(2,at=seq(0,1,1/max((ncol(mm2)-1),1)), label=sub(sample_column_postfix, "", colnames(mm2)), las=2, cex.axis=1, tick=F, family=fontfamily)
    }, mm)
    par(xpd=FALSE);
    legend("bottom", inset=c(0,-.1), legend=c("Homozygous deletion", "Loss", "Gain", "Amplification"),
           fill=c("blue", "lightblue", "darksalmon", "red"), xpd=TRUE, ncol=2)
    par(xpd=TRUE);
    dev.off()
}

geneCN_tab <- read.table(geneCN, sep="\t", header=T, stringsAsFactors=F)
if (!opt$includeChrY) {
    geneCN_tab <- geneCN_tab[geneCN_tab$chrom != "Y",]
}
if (!is.null(opt$summaryConfig)) {
    library(yaml)
    summaryConfig <- yaml.load_file("summary_config.yaml")
    # replace sample names if specified in config
    if ("sample_rename" %in% names(summaryConfig)) {
        for (c in names(summaryConfig$sample_rename)) {
            colnames(geneCN_tab) <- sub(c, summaryConfig$sample_rename[c], colnames(geneCN_tab));
        }
    }
} else {
    summaryConfig <- NULL;
}
plot_heatmap(geneCN_tab, outFile, opt$sampleColumnPostFix, fontfamily, summaryConfig=summaryConfig)
