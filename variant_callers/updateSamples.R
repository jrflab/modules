#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--patient", default = NA, type = 'character', help = "patient id"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

vars = read.csv(file=paste0("sufam/", opt$patient, ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
index = grep("MAF", colnames(vars))
sample_names = unlist(lapply(strsplit(colnames(vars)[index], "_"), function(x) {return(x[2])}))
sample_names[1] = opt$patient

## run sufam
chr = vars$Chromosome
pos = vars$Position
id = rep(".", nrow(vars))
ref = vars$Ref
alt = vars$Alt
qual = rep(100, nrow(vars))
filter = rep("PASS", nrow(vars))
info = rep(".", nrow(vars))
vcf = cbind(chr, pos, id, ref, alt, qual, filter, info)
colnames(vcf) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
write.table(vcf, file=paste0("sufam/", opt$patient, ".vcf"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

for (i in 1:length(sample_names)) {
	system(paste0("source ~/share/usr/anaconda-envs/jrflab-modules-0.1.4/bin/activate ~/share/usr/anaconda-envs/sufam-dev && sufam ~share/reference/GATK_bundle/2.3/human_g1k_v37.fa sufam/", opt$patient, ".vcf bam/", sample_names[i], ".bam"))
}

## fix depth
write.table(vars, file=paste0("sufam/", opt$patient, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
