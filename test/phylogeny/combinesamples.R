#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
				  	make_option("--normal_samples", default = NA, type = 'character', help = "normal samples"),
				  	make_option("--type", default = NA, type = 'character', help = "allele specific or total copy")
				 )
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_samples = na.omit(unlist(strsplit(opt$sample_set, split="_", fixed=TRUE)))
normal_samples = na.omit(unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE)))
normal_samples = normal_samples[normal_samples %in% all_samples]
tumor_samples = all_samples[!(all_samples %in% normal_samples)]

if (opt$type=="allele_specific") {

	CN = list()
	for (i in 1:length(tumor_samples)) {
		load(paste0("facets/cncf/", tumor_samples[i], "_", normal_samples, ".Rdata"))
		CN[[i]] = out2$jointseg[,c("chrom", "maploc", "cnlr", "vafT", "het"),drop=FALSE]
		colnames(CN[[i]]) = c("Chromosome", "Position", "Log2Ratio", "BAF", "Genotype")
	}
	index = lapply(CN, function(x) {paste0(x[,1], ":", x[,2])})
	featureNames = unique(unlist(index))
	for (i in 1:length(index)) {
		featureNames = intersect(featureNames, index[[i]])
	}
	chr = as.numeric(unlist(lapply(strsplit(featureNames, ":", fixed=TRUE), function(x) { x[1] })))
	pos = as.numeric(unlist(lapply(strsplit(featureNames, ":", fixed=TRUE), function(x) { x[2] })))
	index = order(pos, decreasing=FALSE)
	chr = chr[index]
	pos = pos[index]
	index = order(chr, decreasing=FALSE)
	chr = chr[index]
	pos = pos[index]
	featureNames = paste0(chr, ":", pos)
	for (i in 1:length(CN)) {
		rownames(CN[[i]]) = paste0(CN[[i]][,1], ":", CN[[i]][,2])
		CN[[i]] = CN[[i]][featureNames,,drop=FALSE]
	}
	Log2Ratio = do.call(cbind, lapply(CN, function(x) { return(x[,"Log2Ratio"]) } ))
	BAF = do.call(cbind, lapply(CN, function(x) { return(x[,"BAF"]) } ))
	Genotype = do.call(cbind, lapply(CN, function(x) { return(x[,"Genotype"]) } ))
	annotation = data.frame(Chromosome=chr,
							Position=pos)
	colnames(Log2Ratio) = colnames(BAF) = tumor_samples
	save(Log2Ratio, BAF, Genotype, annotation, file=paste0("medicc/allele_specific/mad/", opt$sample_set, ".RData"))
	
} else if (opt$type=="total_copy") {

	CN = list()
	for (i in 1:length(tumor_samples)) {
		load(paste0("facets/cncf/", tumor_samples[i], "_", normal_samples, ".Rdata"))
		CN[[i]] = out2$jointseg[,c("chrom", "maploc", "cnlr"),drop=FALSE]
		colnames(CN[[i]]) = c("Chromosome", "Position", "Log2Ratio")
	}
	index = lapply(CN, function(x) {paste0(x[,1], ":", x[,2])})
	featureNames = unique(unlist(index))
	for (i in 1:length(index)) {
		featureNames = intersect(featureNames, index[[i]])
	}
	chr = as.numeric(unlist(lapply(strsplit(featureNames, ":", fixed=TRUE), function(x) { x[1] })))
	pos = as.numeric(unlist(lapply(strsplit(featureNames, ":", fixed=TRUE), function(x) { x[2] })))
	index = order(pos, decreasing=FALSE)
	chr = chr[index]
	pos = pos[index]
	index = order(chr, decreasing=FALSE)
	chr = chr[index]
	pos = pos[index]
	featureNames = paste0(chr, ":", pos)
	for (i in 1:length(CN)) {
		rownames(CN[[i]]) = paste0(CN[[i]][,1], ":", CN[[i]][,2])
		CN[[i]] = CN[[i]][featureNames,,drop=FALSE]
	}
	Log2Ratio = do.call(cbind, lapply(CN, function(x) { return(x[,"Log2Ratio"]) } ))
	annotation = data.frame(Chromosome=chr,
							Position=pos)
	colnames(Log2Ratio) = tumor_samples
	save(Log2Ratio, annotation, file=paste0("medicc/total_copy/mad/", opt$sample_set, ".RData"))
	
}