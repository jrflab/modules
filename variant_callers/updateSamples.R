#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("CNtu"))
suppressPackageStartupMessages(library("readr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--patient", default = NA, type = 'character', help = "patient id"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

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
vars = cbind(vars, q_t, q_2)
 
#====================================
# loh
#====================================
for (i in 2:length(sample_names)) {
	loh = rep(0, nrow(vars))
	for (j in 1:nrow(vars)) {
		if (q_t[j,i]==q_2[j,i]) {
			loh[j] = 1
		}
	}
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
	load(paste0("ascat/ascat/", sample_names[i], "_", sample_names[1], ".RData"))
	f_hat = vars[,paste0("MAF_", sample_names[i])]
	n = vars[,paste0("DP_", sample_names[i])]
	qt = vars[,paste0("qt_", sample_names[i])]
	qt[qt>10] = 10
	q2 = vars[,paste0("q2_", sample_names[i])]
	q2[q2>10] = 10
	alpha = seq(.1, to=.9, length=50)
	alpha_hat = list()
	indx = f_hat>.2
	if (sum(indx)==0) {
		indx = f_hat>.1
	}
	if (sum(indx)<5) {
		indx = f_hat>0
	}
	for (j in 1:length(alpha)) {
		alpha_hat[[j]] = cancercellFraction(f_hat[indx], n[indx], qt[indx], q2[indx], alpha[j], e=0.01)
	}
	LL = unlist(lapply(alpha_hat, function(x) {sum(x[,"LL"])}))
	
	pdf(file=paste0("sufam/", sample_names[i], "_", sample_names[1], ".pdf"))
	plot(alpha, LL, type="o", col="steelblue", axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
	axis(1, at = NULL, cex.axis = 1.5, padj = 0.25)
    axis(2, at = NULL, cex.axis = 1.5, las = 1)
    mtext(side = 1, text = expression(alpha), line = 4, cex = 1.5)
    mtext(side = 2, text = expression(Sigma~"LL"), line = 4, cex = 1.5)
    index = which.max(LL)
    title(main = paste0("alpha* = ", signif(alpha[index], 3)), cex.main = 1.5)
    box(lwd = 2)
	dev.off()
	
	index = which.max(LL)
	alpha_hat = cancercellFraction(f_hat, n, qt, q2, ifelse((alpha[index]-.25)<=0, alpha[index], alpha[index]-.25), e=0.01)
	cancer_cell_fraction = cbind(cancer_cell_fraction, alpha_hat[,"cancer_cell_frac"])
	ccf_95CI_low = cbind(ccf_95CI_low, alpha_hat[,"ccf_95CI_low"])
	ccf_95CI_high = cbind(ccf_95CI_high, alpha_hat[,"ccf_95CI_high"])
	pr_somatic_clonal = cbind(pr_somatic_clonal, alpha_hat[,"Pr_somatic_clonal"])
	ll = cbind(ll, alpha_hat[,"LL"])
	sq = cbind(sq, alpha_hat[,"sq"])
	clonal_estimate = rep("Subclonal", nrow(vars))
	clonal_estimate[cancer_cell_fraction>.75 | pr_somatic_clonal>.5 | ccf_95CI_low>.9] = "Clonal"
	clonal_status = cbind(clonal_status, clonal_estimate)
}
colnames(cancer_cell_fraction) = colnames(ccf_95CI_low) = colnames(ccf_95CI_high) = colnames(pr_somatic_clonal) = colnames(ll) = colnames(sq) = colnames(clonal_status) = sample_names[2:length(sample_names)]
colnames(cancer_cell_fraction) = paste0("CCF_", colnames(cancer_cell_fraction))
colnames(ccf_95CI_low) = paste0("CCF_95CI_Low_", colnames(ccf_95CI_low))
colnames(ccf_95CI_high) = paste0("CCF_95CI_High_", colnames(ccf_95CI_high))
colnames(pr_somatic_clonal) = paste0("Pr_Somatic_Clonal_", colnames(pr_somatic_clonal))
colnames(ll) = paste0("LL_", colnames(ll))
colnames(sq) = paste0("sq_", colnames(sq))
colnames(clonal_status) = paste0("Clonal_Status_", colnames(clonal_status))

vars = cbind(vars, cancer_cell_fraction,
				   ccf_95CI_low,
				   ccf_95CI_high,
				   pr_somatic_clonal,
				   ll,
				   sq,
				   clonal_status)

write.table(vars, file=paste0("sufam/", opt$patient, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
