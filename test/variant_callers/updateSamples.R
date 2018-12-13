#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("CNtu"))
suppressPackageStartupMessages(library("readr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--patient", default = NA, type = 'character', help = "patient name"))
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = TRUE)
opt = arguments$options

vars = read_tsv(file=paste0("sufam/", opt$patient, ".txt"))
col_names = colnames(vars)
vars = as.data.frame(vars)
colnames(vars) = col_names
index = grep("MAF", colnames(vars))
sample_names = unlist(lapply(strsplit(colnames(vars)[index], "_"), function(x) {return(x[2])}))
sample_names[1] = opt$patient

#====================================
# sufam
#====================================
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
 
#====================================
# dp and maf
#====================================
for (i in 1:length(sample_names)) {
 	if (!file.exists(paste0("sufam/", sample_names[i], ".mat"))) {
 		system(paste0("source ~/share/usr/anaconda/bin/activate ~/share/usr/anaconda-envs/sufam-dev && sufam ~/share/reference/GATK_bundle/2.3/human_g1k_v37.fa sufam/", opt$patient, ".vcf bam/", sample_names[i], ".bam > sufam/", sample_names[i], ".mat"))
 	}
 	tmp = read.csv(file=paste0("sufam/", sample_names[i], ".mat"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
 	## fix depth
 	index = grep("DP", colnames(vars))
 	vars[,index[i]] = tmp[,"cov"]
 	## fix maf
 	index = grep("MAF", colnames(vars))
 	vars[,index[i]] = tmp[,"val_maf"]
}
 
#====================================
# qt and q2
#====================================
q_t = rep(NA, nrow(vars))
q_2 = rep(NA, nrow(vars))
for (i in 2:length(sample_names)) {
	q_t = cbind(q_t, rep(NA, nrow(vars)))
 	q_2 = cbind(q_2, rep(NA, nrow(vars)))
}
colnames(q_t) = colnames(q_2) = c("N", sample_names[2:length(sample_names)])
colnames(q_t) = paste0("qt_", colnames(q_t))
colnames(q_2) = paste0("q2_", colnames(q_2))
vars = cbind(vars, q_t, q_2)
 
#====================================
# loh
#====================================
for (i in 2:length(sample_names)) {
	loh = rep(NA, nrow(vars))
	vars[,paste0("LOH_", sample_names[i])] = loh
}

#====================================
# ccf
#====================================
cancer_cell_fraction = NULL
ccf_95CI_low = NULL
ccf_95CI_high = NULL
pr_somatic_clonal = NULL
ll = NULL
sq = NULL
clonal_status = NULL
for (i in 2:length(sample_names)) {
	cancer_cell_fraction = cbind(cancer_cell_fraction, rep(NA, nrow(vars)))
	ccf_95CI_low = cbind(ccf_95CI_low, rep(NA, nrow(vars)))
	ccf_95CI_high = cbind(ccf_95CI_high, rep(NA, nrow(vars)))
	pr_somatic_clonal = cbind(pr_somatic_clonal, rep(NA, nrow(vars)))
	ll = cbind(ll, rep(NA, nrow(vars)))
	sq = cbind(sq, rep(NA, nrow(vars)))
	clonal_status = cbind(clonal_status, rep(NA, nrow(vars)))
}
colnames(cancer_cell_fraction) = colnames(ccf_95CI_low) = colnames(ccf_95CI_high) = colnames(pr_somatic_clonal) = colnames(ll) = colnames(sq) = colnames(clonal_status) = sample_names[2:length(sample_names)]
colnames(cancer_cell_fraction) = paste0("CCF_", colnames(cancer_cell_fraction))
colnames(ccf_95CI_low) = paste0("CCF_95CI_Low_", colnames(ccf_95CI_low))
colnames(ccf_95CI_high) = paste0("CCF_95CI_High_", colnames(ccf_95CI_high))
colnames(pr_somatic_clonal) = paste0("Pr_Somatic_Clonal_", colnames(pr_somatic_clonal))
colnames(ll) = paste0("LL_", colnames(ll))
colnames(sq) = paste0("sq_", colnames(sq))
colnames(clonal_status) = paste0("Clonal_Status_", colnames(clonal_status))

vars = cbind(vars,
			 cancer_cell_fraction,
			 ccf_95CI_low,
			 ccf_95CI_high,
			 pr_somatic_clonal,
			 ll,
			 sq,
			 clonal_status)

write.table(vars, file=paste0("sufam/", opt$patient, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
