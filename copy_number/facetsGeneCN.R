#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
suppressPackageStartupMessages(library("Cairo"));

optList <- list(
                make_option("--geneLocFile", default = '~/share/reference/IMPACT410_genes_for_copynumber.txt', type = 'character', help = "file containing gene locations"),
                make_option("--outFile", default = NULL, help = "output file"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need facets output files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    facets_files <- arguments$args
}

genes <- read.delim("/home/ngk1/share/reference/IMPACT410_genes_for_copynumber.txt", as.is=T)

genesGR <- GRanges(seqnames=genes$chromosome, 
        ranges=IRanges(as.numeric(genes$start_position), as.numeric(genes$end_position)),
        mcols=genes[,c("order", "Cyt", "hgnc_symbol")])


mm <- do.call("cbind", lapply(facets_files, function(f) {
    tab <- read.delim(f, as.is=T)
    tab$chrom[which(tab$chrom==23)] <- "X"

    tabGR <- GRanges(seqnames=tab$chrom, 
        ranges=IRanges(as.numeric(tab$loc.start), as.numeric(tab$loc.end)),
        mcols=tab[,-c(1:4)])

    fo <- findOverlaps(tabGR, genesGR)
    rr <- ranges(fo, ranges(tabGR), ranges(genesGR))
    df <- cbind(as.data.frame(fo), as.data.frame(rr))

    df <- cbind(df, mcols(genesGR)[df$subjectHits,], mcols(tabGR)[df$queryHits,])

#when genes span multiple segments
    oo <- tapply(df$mcols.cnlr.median, df$subjectHits, function(x){which.max(abs(x))})
    oo <- oo[match(1:409, names(oo))]
    oo[which(is.na(oo))] <- 1

    df <- df[unlist(lapply(1:409, function(x) { which(df$mcols.order==x)[oo[which(names(oo)==x)]]})),]

    ploidy <- table(df$mcols.tcn)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    df$GL <- 0
    df$GL[which(df$mcols.tcn<ploidy)] <- -1
    df$GL[which(df$mcols.tcn==0)] <- -2
    df$GL[which(df$mcols.tcn>ploidy)] <- 1
    df$GL[which(df$mcols.tcn>=ploidy+4)] <- 2

    df <- df[match(genes$order, df$mcols.order),]
    df$GL
}))
colnames(mm) <- facets_files
mm <- cbind(genes, mm)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)


