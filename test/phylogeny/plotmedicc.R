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

init_tree = read.tree(file=paste0("medicc/medicc/", opt$sample_set, "/tree_final.new"))
init_tree = root(init_tree, outgroup="diploid")
boot_tree = foreach (i=1:100) %dopar% {
	if (file.exists(paste0("medicc/boot/", opt$sample_set, "/", str_pad(string=i, width=3, side="left", pad="0"), "/tree_final.new"))) {
		tree = read.tree(file=paste0("medicc/boot/", opt$sample_set, "/", str_pad(string=i, width=3, side="left", pad="0") ,"/tree_final.new"))
		tree = root(tree, outgroup="diploid")
	}
	return(tree)
}
class(init_tree) = "phylo"
class(boot_tree) = "multiPhylo"
node_labels = prop.clades(init_tree, boot_tree, rooted=TRUE)
init_tree$node.label = node_labels
init_tree$tip.label = gsub("_", "-", init_tree$tip.label)

pdf(file=paste0("medicc/plots/", opt$sample_set, ".pdf"), height=7, width=7)
plotTree(tree=init_tree, color="#8CC63F", lwd=3, offset=1)
nodelabels(node=1:init_tree$Nnode+Ntip(init_tree),
		   pie = cbind(as.numeric(init_tree$node.label),100-as.numeric(init_tree$node.label)),
		   piecol = c("goldenrod3","grey85"),
		   cex = 1)
edgelabels(text=paste0(init_tree$edge.length, " "), cex=.75)
dev.off()
