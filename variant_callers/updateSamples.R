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
	if (!file.exists(paste0("sufam/", sample_names[i], ".tmp"))) {
		system(paste0("source ~/share/usr/anaconda/bin/activate ~/share/usr/anaconda-envs/sufam-dev && sufam ~/share/reference/GATK_bundle/2.3/human_g1k_v37.fa sufam/", opt$patient, ".vcf bam/", sample_names[i], ".bam > sufam/", sample_names[i], ".tmp"))
	}
	tmp = read.csv(file=paste0("sufam/", sample_names[i], ".tmp"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
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
q_t = rep(2, nrow(vars))
q_2 = rep(1, nrow(vars))
for (i in 2:length(sample_names)) {
	load(paste0("ascat/ascat/", sample_names[i], "_", sample_names[1], ".RData"))
	Chromosomes = tmp2$SNPpos[tmp3$seg[,"start"],1]
	Chromosomes[Chromosomes==23] = "X"
	Start = tmp2$SNPpos[tmp3$seg[,"start"],2]
	End = tmp2$SNPpos[tmp3$seg[,"end"],2]
	qt = tmp3$seg[,"nA"] + tmp3$seg[,"nB"]
	q2 = apply(tmp3$seg[,c("nA","nB")], 1, max)
	index = rep(NA, nrow(vars))
	for (j in 1:nrow(vars)) {
		indx = which(Chromosomes==vars[j,"Chromosome"] & Start<=vars[j,"Position"] & End>=vars[j,"Position"])
		if (length(indx)!=0) {
			index[j] = indx
		} else {
			index[j] = NA
		}
	}
	q_t = cbind(q_t, qt[index])
	q_2 = cbind(q_2, q2[index])
}
q_t[is.na(q_t)] = 2
q_2[is.na(q_2)] = 1
colnames(q_t) = colnames(q_2) = c("N", sample_names[2:length(sample_names)])
colnames(q_t) = paste0("qt_", colnames(q_t))
colnames(q_2) = paste0("q2_", colnames(q_2))

#====================================
# loh
#====================================
for (i in 2:length(sample_names)) {
	loh = rep(0, nrow(vars))
	for (j in 1:nrow(vars)) {
		if (qt[j,i]==q2[j,i]) {
			loh[j] = 1
		}
	}
	vars[,paste0("LOH_", sample_names[i])] = loh
}

#====================================
# ccf
#====================================


write.table(vars, file=paste0("sufam/", opt$patient, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
