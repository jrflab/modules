#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("phangorn"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
					make_option("--normal_samples", default = NA, type = 'character', help = "normal samples")
				 )

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_samples = na.omit(unlist(strsplit(opt$sample_set, split="_", fixed=TRUE)))
normal_samples = na.omit(unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE)))
normal_samples = normal_samples[normal_samples %in% all_samples]
tumor_samples = all_samples[!(all_samples %in% normal_samples)]

mutation_summary = read_tsv(file=paste0("sufam/", opt$sample_set, ".tsv"), col_types = cols(.default = col_character()))  %>%
 				   type_convert()

mutation_binary = as.data.frame(mutation_summary[,paste0("CALL_", c(tumor_samples, normal_samples)),drop=FALSE])
colnames(mutation_binary) = gsub("CALL_", "", colnames(mutation_binary))

phy_data = as.phyDat(mutation_binary, type="USER", levels=c(0,1))
phy_tree = pratchet(data=phy_data)
phy_tree_w_bl = acctran(tree=phy_tree, data=phy_data)
phy_tree_w_bl = root(phy_tree_w_bl, outgroup=normal_samples)

'bootstrap_data' <- function(x, N=100)
{
	y = list()
	for (i in 1:N) {
		index = sample(1:nrow(x), size=nrow(x), replace=TRUE)
		y[[i]] = x[index,,drop=FALSE]
	}
	return(y)
}


phy_tree_w_bl_boot = list()
mutation_binary_boot = bootstrap_data(x=mutation_binary)
for (i in 1:length(mutation_binary_boot)) {
	phy_data = as.phyDat(mutation_binary_boot[[i]], type="USER", levels=c(0,1))
	phy_tree = pratchet(data=phy_data)
	phy_tree_w_bl_boot[[i]] = acctran(tree=phy_tree, data=phy_data)
	phy_tree_w_bl_boot[[i]] = root(phy_tree_w_bl_boot[[i]], outgroup=normal_samples)
}

class(phy_tree_w_bl) = "phylo"
class(phy_tree_w_bl_boot) = "multiPhylo"
node_labels = prop.clades(phy_tree_w_bl, phy_tree_w_bl_boot, rooted=TRUE)
phy_tree_w_bl$node.label = node_labels
save(list=ls(all=TRUE), file=paste0("pratchet/", opt$sample_set, "/tree_final.RData"))
