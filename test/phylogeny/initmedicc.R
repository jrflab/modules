#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

	load(paste0("medicc/aspcf/", opt$sample_set, ".RData"))
	q1 = qt-q2
	index = !apply(q2, 1, function(x) { any(is.na(x)) }) & !apply(q1, 1, function(x) { any(is.na(x)) })
	q2 = q2[index,,drop=FALSE]
	q1 = q1[index,,drop=FALSE]
	tmp = tmp[index,,drop=FALSE]
	q2[q2>4] = 4
	q1[q1>4] = 4
	desc = cbind(paste0("chrom", unique(tmp[,"Chromosome"])),
				 paste0("major_chr", unique(tmp[,"Chromosome"]), ".fasta"),
				 paste0("minor_chr", unique(tmp[,"Chromosome"]), ".fasta"))
	write.table(desc, file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/desc.txt"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
	for (i in unique(tmp[,"Chromosome"])) {
		cat(">diploid\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/major_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q2)) {
			cat(paste0(">", gsub("-", "_", colnames(q2)[j]), "\n"), file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q2[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		}
	
	
		cat(">diploid\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/minor_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q1)) {
			cat(paste0(">", gsub("-", "_", colnames(q1)[j]), "\n"), file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q1[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/medicc/allele_specific/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		}
	}
	
	load(paste0("medicc/mpcf/", opt$sample_set, ".RData"))
	q2 = round(qt/2)
	q1 = qt - q2
	index = !apply(q2, 1, function(x) { any(is.na(x)) }) & !apply(q1, 1, function(x) { any(is.na(x)) })
	q2 = q2[index,,drop=FALSE]
	q1 = q1[index,,drop=FALSE]
	tmp = tmp[index,,drop=FALSE]
	q2[q2>4] = 4
	q1[q1>4] = 4
	desc = cbind(paste0("chrom", unique(tmp[,"Chromosome"])),
				 paste0("major_chr", unique(tmp[,"Chromosome"]), ".fasta"),
				 paste0("minor_chr", unique(tmp[,"Chromosome"]), ".fasta"))
	write.table(desc, file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/desc.txt"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
	for (i in unique(tmp[,"Chromosome"])) {
		cat(">diploid\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/major_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q2)) {
			cat(paste0(">", gsub("-", "_", colnames(q2)[j]), "\n"), file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q2[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		}
	
	
		cat(">diploid\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/minor_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q1)) {
			cat(paste0(">", gsub("-", "_", colnames(q1)[j]), "\n"), file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q1[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/medicc/total_copy/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		}
	}

