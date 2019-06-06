#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (!dir.exists(paste0("medicc/boot/allele_specific/", opt$sample_set))) {
	dir.create(paste0("medicc/boot/allele_specific/", opt$sample_set))
}

load(paste0("medicc/aspcf/", opt$sample_set, ".RData"))
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

set.seed(0)
for (ii in 1:100) {
	n = nchar(ii)
	if (n==1) {
		n = paste0("00", ii)
	} else if (n==2) {
		n = paste0("0", ii)
	} else {
		n = ii
	}
	index = order(sample(x=1:nrow(tmp), size=nrow(tmp), replace=TRUE))
	q2_b = q2[index,,drop=FALSE]
	q1_b = q1[index,,drop=FALSE]
	tmp_b = tmp[index,,drop=FALSE]
	desc = cbind(paste0("chrom", unique(tmp_b[,"Chromosome"])),
			 	 paste0("major_chr", unique(tmp_b[,"Chromosome"]), ".fasta"),
			 	 paste0("minor_chr", unique(tmp_b[,"Chromosome"]), ".fasta"))
	if (!dir.exists(paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n))) {
		dir.create(paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n))
	}
	write.table(desc, file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/desc.txt"), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
	for (i in unique(tmp[,"Chromosome"])) {
		cat(">diploid\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/major_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/major_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q2_b)) {
			cat(paste0(">", gsub("-", "_", colnames(q2_b)[j]), "\n"), file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/major_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q2_b[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/major_chr", i, ".fasta"), append=TRUE)
		}
		cat(">diploid\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/minor_chr", i, ".fasta"), append=FALSE)
		cat(paste0(rep(1, sum(tmp[,"Chromosome"]==i)), collapse=""), "\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/minor_chr", i, ".fasta"), append=TRUE)
		for (j in 1:ncol(q1_b)) {
			cat(paste0(">", gsub("-", "_", colnames(q1_b)[j]), "\n"), file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/minor_chr", i, ".fasta"), append=TRUE)
			cat(paste0(q1_b[tmp[,"Chromosome"]==i,j], collapse=""), "\n", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/", n, "/minor_chr", i, ".fasta"), append=TRUE)
		}
	}
	if (ii==100) {
		cat("done!", file=paste0("medicc/boot/allele_specific/", opt$sample_set, "/init.timestamp"))
	}
}
