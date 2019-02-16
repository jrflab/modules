#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("RMySQL"))

optList <- list(
				make_option("--sample_names", default = NULL, help = "list of sample names")
			   )

parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options


sample_names = unlist(strsplit(opt$sample_names, split=" ", fixed=TRUE))
genes = read.csv(file="~/share/reference/annotation_gene_lists/geneCN.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes = genes %>%
		filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
		filter(!duplicated(hgnc)) %>%
		arrange(as.integer(chrom), start, end)

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
			
mm <- lapply(facetsFiles, function(f) {
    load(f)
    tab <- fit$cncf
	tab$chrom[which(tab$chrom==23)] <- "X"
	tab$chrom[which(tab$chrom==24)] <- "Y"

	tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
	mcols(tabGR) <- tab %>% select(cnlr.median:lcn.em)

	fo <- findOverlaps(tabGR, genesGR)

	df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

	ploidy <- table(df$tcn.em)
	ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

	df$GL <- 0
	df$GL[df$tcn.em < ploidy] <- -1
	df$GL[df$tcn.em == 0] <- -2
	df$GL[df$tcn.em > ploidy] <- 1
	df$GL[df$tcn.em >= ploidy + 4] <- 2

	load(f)
	noise <- median(abs(out2$jointseg$cnlr-  unlist(apply(out2$out[,c("cnlr.median", "num.mark")], 1, function(x) {rep(x[1], each=x[2])}))))

	lrr <- sort(out2$jointseg$cnlr)
	if (noise <= 0.2) { lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]
	} else if ( noise <= 0.3 ) { lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
	} else { lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]}

	df$GL2 <- 0
	df$GL2[df$cnlr.median < median(lrr)-(2.5*sd(lrr))] <- -1
	df$GL2[df$cnlr.median < median(lrr)-(7*sd(lrr))] <- -2
	df$GL2[df$cnlr.median > median(lrr)+(2*sd(lrr))] <- 1
	df$GL2[df$cnlr.median > median(lrr)+(6*sd(lrr))] <- 2

	df %>% select(hgnc, GL, GL2) %>% ungroup
})
names(mm) <- facetsFiles
for (f in facetsFiles) {
	n <- sub('\\..*', '', sub('.*/', '', f))
	colnames(mm[[f]])[2:3] <- paste(n, c("EM", "LRR_threshold"), sep="_")
}

mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

