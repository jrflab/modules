#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("phytools"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_set", default = NA, type = 'character', help = "sample names set")
				  )

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load(paste0("pratchet/", opt$sample_set, "/tree_final.RData"))

pdf(file=paste0("pratchet/", opt$sample_set, "/tree_final.pdf"), height=7, width=7)
plot.phylo(x=phy_tree_w_bl, edge.color="#8CC63F", edge.width=3, label.offset=1)
nodelabels(node=1:phy_tree_w_bl$Nnode+Ntip(phy_tree_w_bl),
		   pie = cbind(as.numeric(phy_tree_w_bl$node.label),100-as.numeric(phy_tree_w_bl$node.label)),
		   piecol = c("goldenrod3","grey85"),
		   cex = 1)
edgelabels(text=paste0(phy_tree_w_bl$edge.length, " "), cex=.75)
dev.off()
