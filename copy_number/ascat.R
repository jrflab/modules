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

	pdf(file=opt$file_out, width=10, height=4.25)
	par(mar=c(5, 5, 4, 2)+.1)
	CN = out2$jointseg[,c("chrom", "maploc", "cnlr"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
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
	col = rep("grey75", nrow(CN))
	plot(CN[,"Position"], CN[,"Log2Ratio"], type="p", pch=".", cex=1.95, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(CN[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0, col="red")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    rect(xleft=1-1e10, xright=max(CN[,"Position"])+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/log2/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
	dev.off()

} else if (opt$type=="bafall") {
	
	pdf(file=opt$file_out, width=10, height=4.25)
	par(mar=c(5, 5, 4, 2)+.1)
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
	col = rep("grey75", nrow(BAF))
	plot(BAF[,"Position"], BAF[,"BAF"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1.125))
	axis(2, at = NULL, labels = NULL, cex.axis = 1, las = 1)
	mtext(side = 2, text = expression("BAF"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(BAF[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0.5, col="red")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    rect(xleft=1-1e10, xright=max(BAF[,"Position"])+1e10, ybottom=1, ytop=1.25, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/bafall/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
	dev.off()

} else if (opt$type=="bafhet") {

	pdf(file=opt$file_out, width=10, height=4.25)
	par(mar=c(5, 5, 4, 2)+.1)
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
	col = rep("grey75", nrow(BAF))
	plot(BAF[,"Position"], BAF[,"BAF"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1.125))
	axis(2, at = NULL, labels = NULL, cex.axis = 1, las = 1)
	mtext(side = 2, text = expression("BAF"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(BAF[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0.5, col="red")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    rect(xleft=1-1e10, xright=max(BAF[,"Position"])+1e10, ybottom=1, ytop=1.25, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/bafhet/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
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
	col = rep("grey75", nrow(CN_and_BAF))
	pdf(file=opt$file_out, width=10, height=4.25*2)
	par(mar=c(5, 5, 4, 2)+.1, mfrow=c(2,1))
	
	plot(CN_and_BAF[,"Position"], CN_and_BAF[,"Log2Ratio"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
	for (j in 1:nrow(tmp)) {
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
	}
	axis(2, at = NULL, labels = NULL, cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(CN_and_BAF[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0, col="red")
	axis(1, at = .5*(start+end), labels=rep(" ", 23), cex.axis = 0.85, las = 1)
    rect(xleft=1-1e10, xright=max(CN_and_BAF[,"Position"])+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/log2nbaf/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
		
	plot(CN_and_BAF[,"Position"], CN_and_BAF[,"BAF"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(0,1.125))
	points(CN_and_BAF[,"Position"], 1-CN_and_BAF[,"BAF"], type="p", pch=".", cex=1, col=col)
	for (j in 1:nrow(tmp)) {
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"BAF"],2), lty=1, lwd=2.75, col="red")
		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(1-tmp[j,"BAF"],2), lty=1, lwd=2.75, col="red")
	}
	axis(2, at = NULL, labels = NULL, cex.axis = 1, las = 1)
	mtext(side = 2, text = expression("BAF"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(CN_and_BAF[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0.5, col="red")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    rect(xleft=1-1e10, xright=max(CN_and_BAF[,"Position"])+1e10, ybottom=1, ytop=1.25, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/log2nbaf/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
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
	col = rep("grey75", nrow(CN_and_BAF))
	pdf(file=opt$file_out, width=10, height=4.2)
	par(mar=c(5, 5, 4, 2)+.1)
	plot(CN[,"pos"], CN[,"Log2Ratio"], type="p", pch=".", cex=1, col=col, axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
	for (j in 1:nrow(tmp)) {
 		lines(x=c(tmp[j,"Start"], tmp[j,"End"]), y=rep(tmp[j,"Log2Ratio"],2), lty=1, lwd=2.75, col="red")
 	}
	axis(2, at = NULL, labels = NULL, cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	for (j in 1:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=1)
	}
	abline(v=max(CN[,"Position"]), col="goldenrod3", lty=3, lwd=1)
	abline(h=0, col="red")
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    load(gsub(".pdf", ".RData", gsub("total", "ascat", opt$file_out)))
	for (k in 1:8) {
		abline(h=(.75*log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col="brown", lty=3)
		mtext(text=k, side=4, line=.5, at=(.75*log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), las=2, cex=.75, col="brown")
	}
    rect(xleft=1-1e10, xright=max(CN_and_BAF[,"Position"])+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = gsub(".pdf", "", gsub("ascat/total/", "", opt$file_out, fixed=TRUE), fixed=TRUE), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
	dev.off()
	
} else if (opt$type=="plot-chr") {

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
	for (ii in 1:23) {
		pdf(file=paste0(opt$file_out, "/chromosome_", ii, ".pdf"))
		par(mar = c(6.1, 6, 4.1, 3))
		zz = split.screen(figs=matrix(c(0,1,.15,1, 0.065,.975,0.1,.4), nrow=2, ncol=4, byrow=TRUE))
		screen(zz[1])
		start = 0
		end = max(as.numeric(CytoBand[CytoBand[,1]==paste0("chr",ii),2]))
		plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
		index = CN[,"Chromosome"]==i
		z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
		z1 = winsorize(data=z0, tau=3.5, k=15)
		z2 = pcf(data=z1, kmin=100, gamma=100)
		tmp = z2[,c("chrom","start.pos","end.pos","mean")]
		colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
		points(z1[,"pos"], z1[,"cnlr"], type="p", pch=".", cex=1.15, col="grey75")
		for (i in 1:nrow(tmp)) {
			points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
		}
		for (i in 1:(nrow(tmp)-1)) {
			points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
		}
		abline(h=0, lwd=1)
		axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
		mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
		opt = list(gamma=1, rho=purity, psi=ploidy)
		for (k in c(1, 3, 5, 8, 12)) {
			abline(h=.09+(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), col="brown", lty=3)
			mtext(text=k, side=4, line=.5, at=.09+(opt$gamma*log2(((opt$rho)*k + (1-opt$rho)*2)/((opt$rho)*opt$psi + (1-opt$rho)*2))), las=2, cex=.75, col="brown")
		}
		box(lwd=2)
		screen(zz[2])
		assembly = read.csv(file="modules/copy_number/hg19_cytoBandIdeo.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
		plotIdeogram(chrom=ii, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
		close.screen(all.screens=TRUE)
		dev.off()
	}
	
}

warnings()
