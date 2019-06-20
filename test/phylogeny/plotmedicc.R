#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("phytools"))

registerDoMC(12)


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
					make_option("--type", default = NA, type = 'character', help = "allele specific or total copy")
				  )

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (opt$type=="allele_specific") {

	phylo_tree = read.tree(file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/tree_final.new"))
	tip_labels = phylo_tree$tip.label
	index = grep("pad00", tip_labels)
	if (length(index)!=0) {
		phylo_tree = drop.tip(phy=phylo_tree, tip=tip_labels[index], trim.internal=TRUE, rooted=FALSE)
	}
	phylo_tree = root(phylo_tree, outgroup="diploid")

	pdf(file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/tree_final.pdf"), height=7, width=7)
	plotTree(tree=phylo_tree, color="#8CC63F", lwd=3, offset=1)
	edgelabels(text=paste0(phylo_tree$edge.length, " "), cex=.75)
	dev.off()
	
} else if (opt$type=="total_copy") {
	
	phylo_tree = read.tree(file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/tree_final.new"))
	tip_labels = phylo_tree$tip.label
	index = grep("pad00", tip_labels)
	if (length(index)!=0) {
		phylo_tree = drop.tip(phy=phylo_tree, tip=tip_labels[index], trim.internal=TRUE, rooted=FALSE)
	}
	phylo_tree = root(phylo_tree, outgroup="diploid")

	pdf(file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/tree_final.pdf"), height=7, width=7)
	plotTree(tree=phylo_tree, color="#8CC63F", lwd=3, offset=1)
	edgelabels(text=paste0(phylo_tree$edge.length, " "), cex=.75)
	dev.off()
	
}
