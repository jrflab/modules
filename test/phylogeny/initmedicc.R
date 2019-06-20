#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

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

	load(paste0("medicc/allele_specific/aspcf/", opt$sample_set, ".RData"))
	q1 = qt-q2
	index = !apply(q2, 1, function(x) { any(is.na(x)) }) & !apply(q1, 1, function(x) { any(is.na(x)) })
	q2 = q2[index,,drop=FALSE]
	q1 = q1[index,,drop=FALSE]
	tmp = tmp[index,,drop=FALSE]
	q2[q2>4] = 4
	q1[q1>4] = 4
	
	if (ncol(q2)<3) {
		q1x = q1
		colnames(q1x) = paste0(colnames(q1), "_pad00")
		q1 = cbind(q1, q1x)
		q2x = q2
		colnames(q2x) = paste0(colnames(q2), "_pad00")
		q2 = cbind(q2, q2x)
	}
	
	desc = cbind(paste0("chrom", unique(tmp[,"Chromosome"])),
				 paste0("major_chr", unique(tmp[,"Chromosome"]), ".fasta"),
				 paste0("minor_chr", unique(tmp[,"Chromosome"]), ".fasta"))
	write.table(desc, file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/desc.txt"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
	for (i in unique(tmp[,"Chromosome"])) {
		cat(">diploid\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q2)) {
			cat(paste0(">", gsub("-", "_", colnames(q2)[j]), "\n"), file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q2[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		}
	
	
		cat(">diploid\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q1)) {
			cat(paste0(">", gsub("-", "_", colnames(q1)[j]), "\n"), file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q1[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/allele_specific/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		}
	}
} else if (opt$type=="total_copy") {
	
	load(paste0("medicc/total_copy/mpcf/", opt$sample_set, ".RData"))
	ploidy = round(apply(((tmp[,"End"]-tmp[,"Start"])*qt)/sum(tmp[,"End"]-tmp[,"Start"]), 2, sum))
	ploidy[ploidy>=4] = 4
	ploidy[ploidy<=2] = 2
	if (length(unique(ploidy))>1) {
		index = which(ploidy==4)
		
		qt_4n = ceiling(apply(qt[,index,drop=FALSE], 1, mean)/2)*2
		qt_4n[qt_4n==0 & apply(qt[,index,drop=FALSE], 1, mean)!=0] = 1
		qt_2n = round(qt_4n/2)
		qt_2n[qt_2n==0 & apply(qt[,index,drop=FALSE], 1, mean)!=0] = 1
		qt = cbind(qt, diploid_ancestor=qt_2n, tetraploid_ancestor=qt_4n)
		
		q2_4n = ceiling(apply(q2[,index,drop=FALSE], 1, mean)/2)*2
		q2_4n[q2_4n==0 & apply(q2[,index,drop=FALSE], 1, mean)!=0] = 1
		q2_2n = round(q2_4n/2)
		q2_2n[q2_2n==0 & apply(q2[,index,drop=FALSE], 1, mean)!=0] = 1	
		q2 = cbind(q2, diploid_ancestor=q2_2n, tetraploid_ancestor=q2_4n)
		
	}
	
	q1 = qt-q2
	index = !apply(q2, 1, function(x) { any(is.na(x)) }) & !apply(q1, 1, function(x) { any(is.na(x)) })
	q2 = q2[index,,drop=FALSE]
	q1 = q1[index,,drop=FALSE]
	tmp = tmp[index,,drop=FALSE]
	
	q2[q2>4] = 4
	q1[q1>4] = 4
	
	if (ncol(q2)<3) {
		q1x = q1
		colnames(q1x) = paste0(colnames(q1), "_pad00")
		q1 = cbind(q1, q1x)
		q2x = q2
		colnames(q2x) = paste0(colnames(q2), "_pad00")
		q2 = cbind(q2, q2x)
	}
	
	desc = cbind(paste0("chrom", unique(tmp[,"Chromosome"])),
				 paste0("major_chr", unique(tmp[,"Chromosome"]), ".fasta"),
				 paste0("minor_chr", unique(tmp[,"Chromosome"]), ".fasta"))
	write.table(desc, file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/desc.txt"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
	for (i in unique(tmp[,"Chromosome"])) {
		cat(">diploid\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q2)) {
			cat(paste0(">", gsub("-", "_", colnames(q2)[j]), "\n"), file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q2[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/major_chr", i, ".fasta"), append=TRUE)
		}
	
	
		cat(">diploid\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q1)) {
			cat(paste0(">", gsub("-", "_", colnames(q1)[j]), "\n"), file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q1[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/total_copy/medicc/", opt$sample_set, "/minor_chr", i, ".fasta"), append=TRUE)
		}
	}
}
