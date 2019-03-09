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

optList <- list(
				make_option("--sample_names", default = NULL, help = "list of sample names")
			   )

parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = unlist(strsplit(opt$sample_names, split=" ", fixed=TRUE))
genes = read.csv(file="~/share/reference/annotation_gene_lists/annotation_impact_468.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
		filter(Chromosome %in% as.character(c(1:22, "X", "Y"))) %>%
		filter(!duplicated(Gene_Symbol)) %>%
		arrange(as.integer(Chromosome), Start, End)

genes_granges = genes %$%
				GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Gene_Symbol = Gene_Symbol)
mm = lapply(1:length(sample_names), function(i, sample_names, genes, genes_granges) {
	cat(i, "of", length(sample_names), "\n")
    load(paste0("cnvkit/called/", sample_names[i], ".RData"))
	tmp[tmp[,"Chromosome"]==23,"Chromosome"] = "X"
	tmp[tmp[,"Chromosome"]==24,"Chromosome"] = "Y"
	tmp_granges = tmp %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End))
	mcols(tmp_granges) = tmp %>% select(Cat5)
	fo = findOverlaps(tmp_granges, genes_granges)
	x = mcols(genes_granges)[subjectHits(fo),]
	y = mcols(tmp_granges)[queryHits(fo),]
	df = data.frame("Gene_Symbol"=x, "Cat5"=y)
	df = df %>%
		 group_by(Gene_Symbol) %>%
		 top_n(1, abs(Cat5))
	z = as.numeric(df$Cat5)
	names(z) = as.character(df$Gene_Symbol)
	z = z[names(z) %in% genes[,1]]
	res = rep(NA, nrow(genes))
	names(res) = genes[,1]
	res[names(z)] = z
	return(res)
}, sample_names, genes, genes_granges)
bygene = do.call(cbind, mm)
colnames(bygene) = sample_names
bygene = cbind(genes, bygene) %>%
	 	 arrange(as.integer(Chromosome), Start, End)

save(bygene, file="cnvkit/summary/bygene.RData")
write.table(bygene, file="cnvkit/summary/bygene.txt", sep="\t", col.names=TRUE, row.names=FALSE, na="", quote=FALSE)
