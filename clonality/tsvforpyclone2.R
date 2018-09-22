#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--tumor_name", default = NULL, help = "tumor name"),
			   make_option("--normal_name", default = NULL, help = "normal name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options
mutation_summary = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
mutation_summary = subset(mutation_summary, TUMOR_SAMPLE==opt$tumor_name & NORMAL_SAMPLE==opt$normal_name)
index = grepl("mutect", mutation_summary[,"variantCaller"]) | grepl("varscan", mutation_summary[,"variantCaller"])
mutation_summary = mutation_summary[index,,drop=FALSE]
index = mutation_summary[,"Variant_Classification"] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation")
mutation_summary = mutation_summary[index,,drop=FALSE]
index = mutation_summary[,"CHROM"] %in% c(1:23, "X")
mutation_summary = mutation_summary[index,,drop=FALSE]

Chromosome = mutation_summary[,"CHROM"]
Chromosome[Chromosome=="X"] = 23
Chromosome = as.numeric(Chromosome)
Position = as.numeric(mutation_summary[,"POS"])
load(paste0("ascat/ascat/", opt$tumor_name, "_", opt$normal_name, ".RData"))
Chromosomes = tmp2$SNPpos[tmp3$seg[,"start"],1]
Start = tmp2$SNPpos[tmp3$seg[,"start"],2]
End = tmp2$SNPpos[tmp3$seg[,"end"],2]
qt = tmp3$seg[,"nA"] + tmp3$seg[,"nB"]
q2 = apply(tmp3$seg[,c("nA","nB")], 1, max)
index = rep(NA, nrow(mutation_summary))
for (j in 1:nrow(mutation_summary)) {
	indx = which(Chromosomes==Chromosome[j] & Start<=Position[j] & End>=Position[j])
	if (length(indx)!=0) {
		index[j] = indx
	} else {
		index[j] = NA
	}
}
q_t = qt[index]
q_2 = q2[index]
q_t[is.na(q_t)] = 2
q_2[is.na(q_2)] = 1
q_t[q_t==0] = q_2[q_t==0]
q_t[q_t==0] = 2

mutation_id = paste0(mutation_summary[,"SYMBOL"], "_", mutation_summary[,"HGVSp_Short"])
ref_counts = round((1-mutation_summary[,"TUMOR_MAF"])*mutation_summary[,"TUMOR_DP"])
var_counts = round(mutation_summary[,"TUMOR_MAF"]*mutation_summary[,"TUMOR_DP"])
normal_cn = rep(2, nrow(mutation_summary))
minor_cn = rep(0, nrow(mutation_summary))
major_cn = q_t

mutation_summary = cbind(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
colnames(mutation_summary) = c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn")
write.table(mutation_summary, file=paste0("pyclone/", opt$tumor_name, "_", opt$normal_name, "/", opt$tumor_name, "_", opt$normal_name, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
