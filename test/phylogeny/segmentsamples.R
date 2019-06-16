#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_set", default = NA, type = 'character', help = "sample names set"),
					make_option("--normal_samples", default = NA, type = 'character', help = "normal samples"),
					make_option("--gamma", default = NA, type = 'character', help = "segmentation parameter gamma"),
					make_option("--nlog2", default = NA, type = 'character', help = "number of clusters in Log2 ratio"),
					make_option("--nbaf", default = NA, type = 'character', help = "number of clusters in BAF"),
					make_option("--type", default = NA, type = 'character', help = "allele specific or total copy")
				 )
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_samples = na.omit(unlist(strsplit(opt$sample_set, split="_", fixed=TRUE)))
normal_samples = na.omit(unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE)))
normal_samples = normal_samples[normal_samples %in% all_samples]
tumor_samples = all_samples[!(all_samples %in% normal_samples)]

if (opt$type=="allele_specific") {

	load(paste0("medicc/allele_specific/mad/", opt$sample_set, ".RData"))
	gamma = ifelse(is.na(as.numeric(opt$gamma)), 50, as.numeric(opt$gamma))
	nlog2 = ifelse(is.na(as.numeric(opt$nlog2)), 10, as.numeric(opt$nlog2))
	nbaf = ifelse(is.na(as.numeric(opt$nbaf)), 15, as.numeric(opt$nbaf))
	index = apply(Genotype, 1, function(x) {sum(x==1)==length(x)})
	Log2Ratio = Log2Ratio[index,,drop=FALSE]
	BAF = BAF[index,,drop=FALSE]
	annotation = annotation[index,,drop=FALSE]
	colnames(Log2Ratio) = paste0("Log2Ratio_", colnames(Log2Ratio))
	colnames(BAF) = paste0("BAF_", colnames(BAF))
	index = BAF>.5
	BAF[index] = 1 - BAF[index]
	CN_and_BAF = cbind(annotation, Log2Ratio, BAF)
	tmp = NULL
	for (i in 1:23) {
		cn_and_baf = subset(CN_and_BAF, CN_and_BAF[,"Chromosome"]==i)
		x = try(multipcf(data=winsorize(data=cn_and_baf, method="mad", tau=2.5, k=15, verbose=FALSE), gamma=gamma, normalize=FALSE, fast=FALSE, verbose=FALSE), silent=TRUE)
		if (!("try-error" %in% is(x))) {
			colnames(x)[1:5] = c("Chromosome", "Arm", "Start", "End", "N")
			tmp = rbind(tmp, x)
		}
	}
	CN_and_BAF = subset(CN_and_BAF, CN_and_BAF[,"Chromosome"] %in% tmp[,"Chromosome"])
	qt = q2 = matrix(NA, nrow=nrow(tmp), ncol=length(tumor_samples))
	colnames(qt) = colnames(q2) = tumor_samples
	for (i in 1:length(tumor_samples)) {
		ascat = new.env()
		load(paste0("ascat/ascat/", tumor_samples[i], "_", normal_samples, ".RData"), envir=ascat)

		'prunesegments.cn' <- function(x, n=10)
		{
			cnm = matrix(NA, nrow=length(x), ncol=length(x))
			for (j in 1:length(x)) {
				cnm[,j] = abs(2^x[j] - 2^x)
			}
			cnt = hclust(as.dist(cnm), "average")
			cnc = cutree(tree=cnt, k=n)
			for (j in unique(cnc)) {
				indx = which(cnc==j)
				if (length(indx)>2) {
					mcl = mean(x[indx])
					scl = sd(x[indx])
					ind = which(x[indx]<(mcl+1.96*scl) & x[indx]>(mcl-1.96*scl))
					x[indx[ind]] = mean(x[indx[ind]])
				} else {
					x[indx] = mean(x[indx])
				}
			}
			return(x)
		}

		'prunesegments.baf' <- function(x, n=10)
		{
			cnm = matrix(NA, nrow=length(x), ncol=length(x))
			for (j in 1:length(x)) {
				cnm[,j] = abs(2^x[j] - 2^x)
			}
			cnt = hclust(as.dist(cnm), "average")
			cnc = cutree(tree=cnt, k=n)
			for (j in unique(cnc)) {
				indx = which(cnc==j)
				if (length(indx)>2) {
					mcl = mean(x[indx])
					scl = sd(x[indx])
					ind = which(x[indx]<(mcl+1.96*scl) & x[indx]>(mcl-1.96*scl))
					x[indx[ind]] = mean(x[indx[ind]])
				} else {
					x[indx] = mean(x[indx])
				}
			}
			return(x)
		}
		tmp[,paste0("Log2Ratio_", tumor_samples[i])] = prunesegments.cn(x=tmp[,paste0("Log2Ratio_", tumor_samples[i])], n=nlog2)
		tmp[,paste0("BAF_", tumor_samples[i])] = prunesegments.baf(x=tmp[,paste0("BAF_", tumor_samples[i])], n=nbaf)
	
		Tumor_LogR = as.numeric(CN_and_BAF[,paste0("Log2Ratio_", tumor_samples[i])])
		Tumor_BAF = as.numeric(CN_and_BAF[,paste0("BAF_", tumor_samples[i])])
		Tumor_LogR_segmented = rep(tmp[,paste0("Log2Ratio_", tumor_samples[i])], times=tmp[,"N"])
		Tumor_BAF_segmented = rep(tmp[,paste0("BAF_", tumor_samples[i])], times=tmp[,"N"])
		SNPpos = CN_and_BAF[,c("Chromosome", "Position"), drop=FALSE]
		names(Tumor_LogR) = names(Tumor_BAF) = names(Tumor_LogR_segmented) = names(Tumor_BAF_segmented) = rownames(SNPpos) = paste0("chr", CN_and_BAF[,"Chromosome"], ":", CN_and_BAF[,"Position"])
		colnames(SNPpos) = c("chrs", "pos")
		ch = list()
		j = 1
		for (j in 1:length(unique(CN_and_BAF[,"Chromosome"]))) {
			index = which(CN_and_BAF[,"Chromosome"]==(unique(CN_and_BAF[,"Chromosome"]))[j])
			ch[[j]] = index
			j = j + 1
		}
		chr = ch
		chrs = unique(CN_and_BAF[,"Chromosome"])
		gender = "2323"
		sexchromosomes = c(23, 24)
		tmp2 = list(Tumor_LogR=Tumor_LogR,
					Tumor_BAF=Tumor_BAF,
					Tumor_LogR_segmented=Tumor_LogR_segmented,
					Tumor_BAF_segmented=Tumor_BAF_segmented,
					SNPpos=SNPpos,
					chromosomes=ch,
					chrnames=chrs,
					gender=gender,
					sexchromosomes=sexchromosomes)
	
		tmp3 = try(runASCAT(lrr=tmp2$Tumor_LogR,
							baf=tmp2$Tumor_BAF,
							lrrsegmented=tmp2$Tumor_LogR_segmented,
							bafsegmented=tmp2$Tumor_BAF_segmented,
							gender=tmp2$gender,
							SNPpos=tmp2$SNPpos,
							chromosomes=tmp2$chromosomes,
							chrnames=tmp2$chrnames,
							sexchromosomes=tmp2$sexchromosomes,
							failedqualitycheck=FALSE,
							distance = paste0("medicc/allele_specific/ascat/", tumor_samples[i], "_", normal_samples, ".pdf"),
							copynumberprofile = NULL,
							nonroundedprofile = NULL, 
							aberrationreliability = NULL,
							gamma = 1, rho_manual = ascat$tmp3$rho, psi_manual = ascat$tmp3$psi, y_limit = 3, circos = NA))
						
		if (!("try-error" %in% is(tmp3))) {
			chr = SNPpos[tmp3$seg_raw[,1],1]
			pos = SNPpos[tmp3$seg_raw[,1],2]
			qt[tmp[,1] %in% chr & tmp[,3] %in% pos,tumor_samples[i]] = tmp3$seg_raw[,"nA"] + tmp3$seg_raw[,"nB"]
			q2[tmp[,1] %in% chr & tmp[,3] %in% pos,tumor_samples[i]] = apply(tmp3$seg_raw[,c("nA", "nB"),drop=FALSE], 1, max, na.rm=TRUE)
		}
	}
	save(list=ls(all=TRUE), file=paste0("medicc/allele_specific/aspcf/", opt$sample_set, ".RData"))
	
} else if (opt$type=="total_copy") {

	load(paste0("medicc/total_copy/mad/", opt$sample_set, ".RData"))
	gamma = ifelse(is.na(as.numeric(opt$gamma)), 150, as.numeric(opt$gamma))
	nlog2 = ifelse(is.na(as.numeric(opt$nlog2)), 10, as.numeric(opt$nlog2))
	colnames(Log2Ratio) = paste0("Log2Ratio_", colnames(Log2Ratio))
	CN_and_BAF = cbind(annotation, Log2Ratio)
	tmp = NULL
	for (i in 1:23) {
		cn_and_baf = subset(CN_and_BAF, CN_and_BAF[,"Chromosome"]==i)
		x = try(multipcf(data=winsorize(data=cn_and_baf, method="mad", tau=2.5, k=15, verbose=FALSE), gamma=gamma, normalize=FALSE, fast=FALSE, verbose=FALSE), silent=TRUE)
		if (!("try-error" %in% is(x))) {
			colnames(x)[1:5] = c("Chromosome", "Arm", "Start", "End", "N")
			tmp = rbind(tmp, x)
		}
	}
	CN_and_BAF = subset(CN_and_BAF, CN_and_BAF[,"Chromosome"] %in% tmp[,"Chromosome"])
	qt = q2 = matrix(NA, nrow=nrow(tmp), ncol=length(tumor_samples))
	colnames(qt) = colnames(q2) = tumor_samples
	for (i in 1:length(tumor_samples)) {
		ascat = new.env()
		load(paste0("ascat/ascat/", tumor_samples[i], "_", normal_samples, ".RData"), envir=ascat)

		'prunesegments.cn' <- function(x, n=10)
		{
			cnm = matrix(NA, nrow=length(x), ncol=length(x))
			for (j in 1:length(x)) {
				cnm[,j] = abs(2^x[j] - 2^x)
			}
			cnt = hclust(as.dist(cnm), "average")
			cnc = cutree(tree=cnt, k=n)
			for (j in unique(cnc)) {
				indx = which(cnc==j)
				if (length(indx)>2) {
					mcl = mean(x[indx])
					scl = sd(x[indx])
					ind = which(x[indx]<(mcl+1.96*scl) & x[indx]>(mcl-1.96*scl))
					x[indx[ind]] = mean(x[indx[ind]])
				} else {
					x[indx] = mean(x[indx])
				}
			}
			return(x)
		}
		
		'absolute.cn' <- function(rho, psi, gamma=1, x)
		{
			rho = ifelse(is.na(rho), 1, rho)
			psi = ifelse(is.na(psi), 2, psi)
			return(invisible(((((2^(x/gamma))*(rho*psi+(1-rho)*2)) - ((1-rho)*2))/rho)))
		}

		tmp[,paste0("Log2Ratio_", tumor_samples[i])] = prunesegments.cn(x=tmp[,paste0("Log2Ratio_", tumor_samples[i])], n=nlog2)
		purity = ifelse(is.na(ascat$tmp3$rho), 1, ascat$tmp3$rho)
		ploidy = ifelse(is.na(ascat$tmp3$psi), 1, ascat$tmp3$psi)
		qt[,tumor_samples[i]] = ifelse(round(absolute.cn(rho=purity, psi=ploidy, x=tmp[,paste0("Log2Ratio_", tumor_samples[i])]))<0, 0, round(absolute.cn(rho=purity, psi=ploidy, x=tmp[,paste0("Log2Ratio_", tumor_samples[i])])))
		q2[,tumor_samples[i]] = ceiling(qt[,tumor_samples[i]]/2)
	}
	save(list=ls(all=TRUE), file=paste0("medicc/total_copy/mpcf/", opt$sample_set, ".RData"))
}
