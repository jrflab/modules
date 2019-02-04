#!/usr/bin/env Rscript

all_vars = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
CHROM = all_vars[,"CHROM"]
POS = all_vars[,"POS"]
ID = all_vars[,"ID"]
REF = all_vars[,"REF"]
ALT = all_vars[,"ALT"]
QUAL = FILTER = rep(".", nrow(all_vars))
INFO = paste0(all_vars[,"SYMBOL"], all_vars[,"HGVSp_Short"])
vcf = data.frame(CHROM, POS, ID, REF, ALT, QUAL, INFO)

cat("#", file="sufam/pdx.vcf", append=FALSE)
write.table(vcf, file="sufam/pdx.vcf", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
