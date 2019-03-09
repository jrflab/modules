#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

vcf_snp = read.csv(file=paste0("vcf_ann/", opt$sample_name, ".gatk_snps.vcf"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
vcf_indel = read.csv(file=paste0("vcf_ann/", opt$sample_name, ".gatk_indels.vcf"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
vcf = rbind(vcf_snp, vcf_indel)
pos = as.numeric(vcf[,2])
index = order(pos, decreasing=FALSE)
vcf = vcf[index,,drop=FALSE]
chr = as.character(vcf[,1])
chr[chr=="X"] = 23
chr[chr=="Y"] = 24
chr = as.numeric(chr)
index = is.na(chr)
chr = chr[!index]
vcf = vcf[!index,,drop=FALSE]
index = order(chr, decreasing=FALSE)
vcf = vcf[index,,drop=FALSE]
vcf = vcf[,1:7,drop=FALSE]
colnames(vcf) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
vcf = cbind(vcf, "INFO"=rep(".", nrow(vcf)))
index = grepl(",", vcf[,"REF"]) | grepl(",", vcf[,"ALT"])
vcf = vcf[!index,,drop=FALSE]
index = duplicated(paste0(vcf[,1], ":", vcf[,2]))
vcf = vcf[!index,,drop=FALSE]

cat("##fileformat=VCFv4.1\n", file=paste0("cravat/", opt$sample_name, ".vcf"), append=FALSE)
write.table(vcf, file=paste0("cravat/", opt$sample_name, ".vcf"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=TRUE)
