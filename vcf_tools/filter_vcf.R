#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

vcf = read.csv(file=paste0("cravat/", opt$sample_name, ".vcf"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
maf = read.csv(file=paste0("cravat/", opt$sample_name, ".maf"), header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
index = maf[,"Variant_Classification"] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")
vcf = vcf[index,,drop=FALSE]
vcf[,1] = paste0("chr", vcf[,1])
colnames(vcf) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
cat("##fileformat=VCFv4.1\n", file=paste0("cravat/", opt$sample_name, ".cravat.vcf"), append=FALSE)
write.table(vcf, file=paste0("cravat/", opt$sample_name, ".cravat.vcf"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=TRUE)
