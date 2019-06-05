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

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

allele_specific_tree = read.tree(file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/tree_final.new"))
allele_specific_tree = root(allele_specific_tree, outgroup="diploid")

pdf(file=paste0("medicc/plots/", opt$sample_set, "_allele_specific.pdf"), height=7, width=7)
plotTree(tree=allele_specific_tree, color="#8CC63F", lwd=3, offset=1)
edgelabels(text=paste0(allele_specific_tree$edge.length, " "), cex=.75)
dev.off()

total_copy_tree = read.tree(file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/tree_final.new"))
total_copy_tree = root(total_copy_tree, outgroup="diploid")

pdf(file=paste0("medicc/plots/", opt$sample_set, "_total_copy.pdf"), height=7, width=7)
plotTree(tree=total_copy_tree, color="#8CC63F", lwd=3, offset=1)
edgelabels(text=paste0(total_copy_tree$edge.length, " "), cex=.75)
dev.off()
