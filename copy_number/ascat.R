#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--type", default = NA, type = 'character', help = "type of analysis"),
				  make_option("--file_in", default = NA, type = 'character', help = "input file name"),
				  make_option("--file_out", default = NA, type = 'character', help = "output file name"),
				  make_option("--gamma", default = NA, type = 'numeric', help = "gamma parameter in pcf"),
				  make_option("--nlog2", default = NA, type = 'numeric', help = "number of clusters in Log2 ratio"),
				  make_option("--nbaf", default = NA, type = 'numeric', help = "number of clusters in BAF"),
				  make_option("--rho", default = NA, type = 'numeric', help = "purity for ASCAT"),
				  make_option("--psi", default = NA, type = 'numeric', help = "ploidy for ASCAT"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options
load("modules/copy_number/CytoBand.RData")
load(opt$file_in)

if (opt$type=="log2") {

	pdf(file=opt$file_out, height=7*10/7*2/3, width=7*20/7)
	CN = out2$jointseg[,c("chrom", "maploc", "cnlr"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
	tmp = winsorize(data=CN[,c("Chromosome","Position","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		CN[CN[,"Chromosome"]==j,"Position"] = CN[CN[,"Chromosome"]==j,"Position"] + start[j]
	}
	col = rainbow_hcl(nrow(CN), start=30, end=300)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(CN[,"Position"], tmp[,"Log2Ratio"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2,2))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(CN[,"Position"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	dev.off()

} else if (opt$type=="bafall") {
	
	pdf(file=opt$file_out, height=7*10/7*2/3, width=7*20/7)
	BAF = out2$jointseg[,c("chrom", "maploc", "vafT"),drop=FALSE]
	colnames(BAF) = c("Chromosome", "Position", "BAF")
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		BAF[BAF[,"Chromosome"]==j,"Position"] = BAF[BAF[,"Chromosome"]==j,"Position"] + start[j]
	}
	col = rainbow_hcl(nrow(BAF), start=30, end=300)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(BAF[,"Position"], BAF[,"BAF"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = "B Allele Frequency", line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0.5, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(BAF[,"Position"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	dev.off()

} else if (opt$type=="bafhet") {

	pdf(file=opt$file_out, height=7*10/7*2/3, width=7*20/7)
	BAF = out2$jointseg[,c("chrom", "maploc", "vafT"),drop=FALSE]
	index = out2$jointseg[,"het"]==1
	BAF = BAF[index,,drop=FALSE]
	colnames(BAF) = c("Chromosome", "Position", "BAF")
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		BAF[BAF[,"Chromosome"]==j,"Position"] = BAF[BAF[,"Chromosome"]==j,"Position"] + start[j]
	}
	col = rainbow_hcl(nrow(BAF), start=30, end=300)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(BAF[,"Position"], BAF[,"BAF"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1))
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = "B Allele Frequency", line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0.5, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(BAF[,"Position"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	dev.off()

} else if (opt$type=="aspcf") {

	gamma = ifelse(is.na(as.numeric(opt$gamma)), 70, as.numeric(opt$gamma))
	
	CN_and_BAF = out2$jointseg[,c("chrom", "maploc", "cnlr", "vafT"),drop=FALSE]
	index = out2$jointseg[,"het"]==1
	CN_and_BAF = CN_and_BAF[index,,drop=FALSE]
	colnames(CN_and_BAF) = c("Chromosome", "Position", "Log2Ratio", "BAF")
	index = CN_and_BAF[,"BAF"]>0.5
	CN_and_BAF[index,"BAF"] = 1 - CN_and_BAF[index,"BAF"]
	tmp = multipcf(data=winsorize(data=CN_and_BAF, method="mad", tau=2.5, k=25, verbose=FALSE), gamma=gamma, fast=FALSE, verbose=FALSE)
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio", "BAF")
	save(CN_and_BAF, tmp, file=opt$file_out)

} else if (opt$type=="plot-aspcf") {

	nlog2 = ifelse(is.na(as.numeric(opt$nlog2)), 10, as.numeric(opt$nlog2))
	nbaf = ifelse(is.na(as.numeric(opt$nbaf)), 15, as.numeric(opt$nbaf))

	'prunesegments.cn' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"Log2Ratio"])
				scl = sd(x[indx,"Log2Ratio"])
				ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
				x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
			} else {
				x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
			}
		}
		return(x)
	}

	'prunesegments.baf' <- function(x, n=15)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"BAF"] - 2^x[,"BAF"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"BAF"])
				scl = sd(x[indx,"BAF"])
				ind = which(x[indx,"BAF"]<(mcl+1.96*scl) & x[indx,"BAF"]>(mcl-1.96*scl))
				x[indx[ind],"BAF"] = mean(x[indx[ind],"BAF"])
			} else {
				x[indx,"BAF"] = mean(x[indx,"BAF"])
			}
		}
		return(x)
	}
	tmp = prunesegments.cn(x=tmp, n=nlog2)
	tmp = prunesegments.baf(x=tmp, n=nbaf)
	
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		CN_and_BAF[CN_and_BAF[,"Chromosome"]==j,"Position"] = CN_and_BAF[CN_and_BAF[,"Chromosome"]==j,"Position"] + start[j]
		tmp[tmp[,"Chromosome"]==j,"Start"] = tmp[tmp[,"Chromosome"]==j,"Start"] + start[j]
		tmp[tmp[,"Chromosome"]==j,"End"] = tmp[tmp[,"Chromosome"]==j,"End"] + start[j]
	}
	col = "grey80"
	pdf(file=opt$file_out, height=2*7*10/7*2/3, width=7*20/7)
	par(mar=c(5, 5, 4, 2)+.1, mfrow=c(2,1))
	plot(CN_and_BAF[,"Position"], CN_and_BAF[,"Log2Ratio"], type="p", pch=".", cex=2.5, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2,2))
	for (j in 1:nrow(tmp)) {
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
	}
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(CN_and_BAF[,"Position"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=rep(" ", 23), cex.axis = 0.85, las = 1)
	plot(CN_and_BAF[,"Position"], CN_and_BAF[,"BAF"], type="p", pch=".", cex=2.5, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1))
	points(CN_and_BAF[,"Position"], 1-CN_and_BAF[,"BAF"], type="p", pch=".", cex=2.5, col=col)
	for (j in 1:nrow(tmp)) {
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"BAF"],2), lty=1, lwd=2.75, col="red")
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(1-tmp[j,"BAF"],2), lty=1, lwd=2.75, col="red")
	}
	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 2, text = "B Allele Frequency", line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0.5, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(CN_and_BAF[,"Position"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	dev.off()
	
} else if (opt$type=="run-ascat") {
	
	nlog2 = ifelse(is.na(as.numeric(opt$nlog2)), 10, as.numeric(opt$nlog2))
	nbaf = ifelse(is.na(as.numeric(opt$nbaf)), 15, as.numeric(opt$nbaf))
	rho = as.numeric(opt$rho)
	psi = as.numeric(opt$psi)
	if (is.na(rho) | is.na(psi)) {
		rho = NA
		psi = NA
	}

	'prunesegments.cn' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"Log2Ratio"])
				scl = sd(x[indx,"Log2Ratio"])
				ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
				x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
			} else {
				x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
			}
		}
		return(x)
	}

	'prunesegments.baf' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"BAF"] - 2^x[,"BAF"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"BAF"])
				scl = sd(x[indx,"BAF"])
				ind = which(x[indx,"BAF"]<(mcl+1.96*scl) & x[indx,"BAF"]>(mcl-1.96*scl))
				x[indx[ind],"BAF"] = mean(x[indx[ind],"BAF"])
			} else {
				x[indx,"BAF"] = mean(x[indx,"BAF"])
			}
		}
		return(x)
	}
	tmp = prunesegments.cn(x=tmp, n=nlog2)
	tmp = prunesegments.baf(x=tmp, n=nbaf)
	
	Tumor_LogR = as.numeric(CN_and_BAF[,"Log2Ratio"])
	Tumor_BAF = as.numeric(CN_and_BAF[,"BAF"])
	Tumor_LogR_segmented = rep(tmp[,"Log2Ratio"], times=tmp[,"N"])
	Tumor_BAF_segmented = rep(tmp[,"BAF"], times=tmp[,"N"])
	SNPpos = CN_and_BAF[,c("Chromosome", "Position"), drop=FALSE]
	names(Tumor_LogR) = names(Tumor_BAF) = names(Tumor_LogR_segmented) = names(Tumor_BAF_segmented) = rownames(SNPpos) = paste0("chr", CN_and_BAF[,"Chromosome"], ":", CN_and_BAF[,"Position"])
	colnames(SNPpos) = c("chrs", "pos")
	ch = list()
	for (j in 1:23) {
		index = which(CN_and_BAF[,"Chromosome"]==j)
		ch[[j]] = index
	}
	chr = ch
	chrs = 1:23
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
        	                distance = opt$file_out,
        	                copynumberprofile = NULL,
        	                nonroundedprofile = NULL, 
        	                aberrationreliability = NULL,
        	                gamma = 1, rho_manual = rho, psi_manual = psi, y_limit = 3, circos = NA))
                        
    if (!("try-error" %in% is(tmp3))) {
        purity = tmp3$rho
        ploidy = sum((tmp3$seg_raw[,"nAraw"]+tmp3$seg_raw[,"nBraw"])*(tmp2$SNPpos[tmp3$seg_raw[,"end"],"pos"]-tmp2$SNPpos[tmp3$seg_raw[,"start"],"pos"])/sum(as.numeric(tmp2$SNPpos[tmp3$seg_raw[,"end"],"pos"]-tmp2$SNPpos[tmp3$seg_raw[,"start"],"pos"])))
    	save(tmp, tmp2, tmp3, CN_and_BAF, purity, ploidy, file=gsub(".pdf", ".RData", opt$file_out))
    }
	
} else if (opt$type=="total-copy") {

	'prunesegments.cn' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"Log2Ratio"])
				scl = sd(x[indx,"Log2Ratio"])
				ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
				x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
			} else {
				x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
			}
		}
		return(x)
	}

	CN = out2$jointseg[,c("chrom", "maploc", "cnlr"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
	tmp = pcf(data=winsorize(data=CN, method="mad", tau=2.5, k=25, verbose=FALSE), kmin = 250, gamma=250, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	save(CN, tmp, file=gsub(".pdf", ".RData", opt$file_out))
	tmp = prunesegments.cn(x=tmp, n=7)
	CN = winsorize(data=CN[,c("Chromosome","Position","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		tmp[tmp[,"Chromosome"]==j,"Start"] = tmp[tmp[,"Chromosome"]==j,"Start"] + start[j]
		tmp[tmp[,"Chromosome"]==j,"End"] = tmp[tmp[,"Chromosome"]==j,"End"] + start[j]
		CN[CN[,"chrom"]==j,"pos"] = CN[CN[,"chrom"]==j,"pos"] + start[j]
	}
	col = "grey80"
	pdf(file=opt$file_out, height=7*10/7*2/3, width=7*20/7)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(CN[,"pos"], CN[,"Log2Ratio"], type="p", pch=".", cex=1.5, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,4))
 	for (j in 1:nrow(tmp)) {
 		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
 	}
 	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 1, text = "Chromosome", line = 3, cex = 1.25)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3")
	abline(h=0, col="red")
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3")
	}
	abline(v=max(CN[,"pos"]), col="goldenrod3")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	
	load(gsub(".pdf", ".RData", gsub("total", "ascat", opt$file_out)))
	
	for (k in 1:8) {
		abline(h=(.75*log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col="darkorange", lty=3)
		mtext(text=k, side=4, line=.5, at=(.75*log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), las=2, cex=.75, col="darkorange")
	}
	dev.off()
	
}

warnings()
