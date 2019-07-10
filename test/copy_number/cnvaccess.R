#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {

	'plot_log2_' <- function(x, title = "") {
		par(mar=c(5, 5, 4, 2)+.1)
		data("CytoBand")
		end = NULL
		for (i in 1:23) {
			end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
		}
		end = cumsum(end)
		start = c(1, end[1:22]+1)
		CytoBand = cbind(start, end)
		index = NULL
		for (i in 1:23) {
			index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chromosome==i)))
		}
		plot(index, x$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		abline(v=1, col="goldenrod3", lty=3, lwd=.5)
		abline(h=0, col="red", lty=1, lwd=1)
		for (j in 1:23) {
			abline(v=CytoBand[j,"end"], col="goldenrod3", lty=3, lwd=.5)
		}
		axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
		rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
		title(main = title, line=-1, cex.main=.75, font.main=1)
		box(lwd=1.5)
	}
	
	dataA = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-A.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataA[dataA[,"chromosome"]=="X", "chromosome"] = 23
	dataA[dataA[,"chromosome"]=="Y", "chromosome"] = 24
	dataA[,"chromosome"] = as.numeric(dataA[,"chromosome"])
	dataA = subset(dataA, dataA[,"chromosome"]<=23)
	
	tmp = dataA[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataA[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataB = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-B.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataB[dataB[,"chromosome"]=="X", "chromosome"] = 23
	dataB[dataB[,"chromosome"]=="Y", "chromosome"] = 24
	dataB[,"chromosome"] = as.numeric(dataB[,"chromosome"])
	dataB = subset(dataB, dataB[,"chromosome"]<=23)
	
	tmp = dataB[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataB[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataC = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".no-pool.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataC[dataC[,"chromosome"]=="X", "chromosome"] = 23
	dataC[dataC[,"chromosome"]=="Y", "chromosome"] = 24
	dataC[,"chromosome"] = as.numeric(dataC[,"chromosome"])
	dataC = subset(dataC, dataC[,"chromosome"]<=23)
	
	tmp = dataC[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=1.5, k=25, verbose=FALSE, return.outliers=TRUE)
	dataC[tmp$wins.outliers[,3]!=0,"log2"] = NA

	data = rbind(dataA, dataB, dataC)
	index = order(data[,2])
	data = data[index,,drop=FALSE]
	index = order(data[,1])
	data = data[index,,drop=FALSE]
	
	pdf(file = paste0("cnvaccess/plot/log2/", opt$sample_name, ".pdf"), width=10, height=4.2)
	plot_log2_(x=data, title=opt$sample_name)
	dev.off()
	
} else if (as.numeric(opt$type)==2) {

	data("CytoBand")

	dataA = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-A.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataA[dataA[,"chromosome"]=="X", "chromosome"] = 23
	dataA[dataA[,"chromosome"]=="Y", "chromosome"] = 24
	dataA[,"chromosome"] = as.numeric(dataA[,"chromosome"])
	dataA = subset(dataA, dataA[,"chromosome"]<=23)
	
	tmp = dataA[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataA[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataB = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-B.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataB[dataB[,"chromosome"]=="X", "chromosome"] = 23
	dataB[dataB[,"chromosome"]=="Y", "chromosome"] = 24
	dataB[,"chromosome"] = as.numeric(dataB[,"chromosome"])
	dataB = subset(dataB, dataB[,"chromosome"]<=23)
	
	tmp = dataB[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataB[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataC = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".no-pool.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataC[dataC[,"chromosome"]=="X", "chromosome"] = 23
	dataC[dataC[,"chromosome"]=="Y", "chromosome"] = 24
	dataC[,"chromosome"] = as.numeric(dataC[,"chromosome"])
	dataC = subset(dataC, dataC[,"chromosome"]<=23)
	
	tmp = dataC[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=1.5, k=25, verbose=FALSE, return.outliers=TRUE)
	dataC[tmp$wins.outliers[,3]!=0,"log2"] = NA

	CN = rbind(dataA, dataB, dataC)
	index = order(CN[,2])
	CN = CN[index,,drop=FALSE]
	index = order(CN[,1])
	CN = CN[index,,drop=FALSE]
	CN = CN[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")

	
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
	
	'plotIdeogram' <- function (chrom, cyto.text = FALSE, cex = 0.6, cyto.data, cyto.unit = "bp", unit) {
		if (chrom == 23) {
        		chrom.cytoband <- cyto.data[cyto.data[, 1] == "chrX", ]
    	} else {
        	if (chrom == 24) {
            	chrom.cytoband <- cyto.data[cyto.data[, 1] == "chrY", ]
        	} else {
            	chrom.cytoband <- cyto.data[cyto.data[, 1] == paste("chr", chrom, sep = ""), ]
        	}
    	}
    	cyto.start <- chrom.cytoband[, 2]
    	cyto.end <- chrom.cytoband[, 3]
    	scale <- copynumber:::convert.unit(unit1 = unit, unit2 = cyto.unit)
    	xleft <- cyto.start * scale
    	xright <- cyto.end * scale
    	n <- length(xleft)
    	chrom.length <- xright[n] - xleft[1]
    	stain <- chrom.cytoband[, 5]
    	sep.stain <- c("gpos", "gneg", "acen", "gvar", "stalk")
    	g <- sapply(sep.stain, grep, x = stain, fixed = TRUE)
    	centromere <- g$acen
    	stalk <- g$stalk
    	col <- rep("", n)
    	col[stain == "gneg"] <- "white"
    	col[stain == "gpos100"] <- "black"
    	col[stain == "gpos75"] <- "gray25"
    	col[stain == "gpos50"] <- "gray50"
    	col[stain == "gpos25"] <- "gray75"
    	col[stain == "stalk"] <- "gray90"
    	col[stain == "gvar"] <- "grey"
    	col[stain == "acen"] <- "yellow"
    	density <- rep(NA, n)
    	angle <- rep(45, n)
    	density[stain == "gvar"] <- 15
    	ylow <- 0
    	yhigh <- 1
    	plot(x = c(0, max(xright)), y = c(ylow, yhigh), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, max(xright)), ylim = c(0, 1), xaxs = "i")
    	skip.rect <- c(1, n, stalk)
    	rect(xleft[-skip.rect], rep(ylow, n - length(skip.rect)), xright[-skip.rect], rep(yhigh, n - length(skip.rect)), 
        col = col[-skip.rect], border = "black", density = density[-skip.rect], 
        angle = angle[-skip.rect])
    	draw.roundEdge(start = xleft[1], stop = xright[1], y0 = ylow, y1 = yhigh, col = col[1], bow = "left", density = density[1], angle = angle[1], chrom.length = chrom.length)
    	draw.roundEdge(start = xleft[n], stop = xright[n], y0 = ylow, y1 = yhigh, col = col[n], bow = "right", density = density[n], angle = angle[n], chrom.length = chrom.length)
    	if (length(stalk) > 0) {
        	for (i in 1:length(stalk)) {
            	copynumber:::drawStalk(xleft[stalk[i]], xright[stalk[i]], ylow, yhigh, col = col[stalk[i]])
        	}
    	}
    	if (cyto.text) {
    		mtext(text = paste(chrom.cytoband[, 4], "-", sep = " "), side = 1, at = (xleft + (xright - xleft)/2), cex = cex, las = 2, adj = 1, xpd = NA)
    	}
	}

	'draw.roundEdge' <- function (start, stop, y0, y1, col, bow, density = NA, angle = 45, lwd = 1, chrom.length) {
    	f <- rep(0, 0)
    	f[1] <- 0.001
    	i = 1
    	half <- y0 + (y1 - y0)/2
    	while (f[i] < half) {
    	    f[i + 1] <- f[i] * 1.3
    	    i <- i + 1
    	}
    	f <- f[-length(f)]
    	Y <- c(y1, y1, y1 - f, half, y0 + rev(f), y0, y0)
    	cyto.length <- stop - start
    	share <- cyto.length/chrom.length
    	if (share > 0.2) {
    	    share <- 0.2
    	}
    	if (bow == "left") {
    	    round.start <- start + cyto.length * (1 - share)^20
    	    x <- seq(round.start, start, length.out = (length(f) + 2))
        	revx <- rev(x[-length(x)])
        	x <- c(x, revx)
        	X <- c(stop, x, stop)
    	} else {
        	if (bow == "right") {
            	round.start <- stop - cyto.length * (1 - share)^20
            	x <- seq(round.start, stop, length.out = (length(f) + 2))
            	revx <- rev(x[-length(x)])
            	x <- c(x, revx)
            	X <- c(start, x, start)
        	}
    	}
    	polygon(x = X, y = Y, col = col, border = "black", density = density, angle = angle, lwd = lwd)
	}

	for (ii in 1:23) {
		
		pdf(file=paste0("cnvaccess/plot/bychr/", opt$sample_name, "/chromosome_", ii, ".pdf"))
		par(mar = c(6.1, 6, 4.1, 3))
		zz = split.screen(figs=matrix(c(0,1,.15,1, 0,1,0.0775,.4), nrow=2, ncol=4, byrow=TRUE))
		screen(zz[1])
		start = 0
		end = max(as.numeric(CytoBand[CytoBand[,1]==ii,4]))
		plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
		index = CN[,"Chromosome"]==ii
		z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
		z1 = pcf(data=z0, kmin=10, gamma=50)
		tmp = z1[,c("chrom","start.pos","end.pos","mean")]
		colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
		points(z0[,"Position"], z0[,"Log2Ratio"], type="p", pch=".", cex=2, col="grey80")
		for (i in 1:nrow(tmp)) {
			points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
		}
		for (i in 1:(nrow(tmp)-1)) {
			points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
		}
		abline(h=0, lwd=1)
		axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
		mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
		# load(paste0(gsub("bychr", "ascat", opt$file_out, fixed=TRUE), ".RData"))
		# z3 = list(gamma=1, rho=purity, psi=ploidy)
		# for (k in c(1, 3, 5, 8, 12)) {
		#	 abline(h=.09+(z3$gamma*log2(((z3$rho)*k + (1-z3$rho)*2)/((z3$rho)*z3$psi + (1-z3$rho)*2))), col="brown", lty=3)
		#	 mtext(text=k, side=4, line=.5, at=.09+(z3$gamma*log2(((z3$rho)*k + (1-z3$rho)*2)/((z3$rho)*z3$psi + (1-z3$rho)*2))), las=2, cex=.75, col="brown")
		# }
		box(lwd=2)
		screen(zz[2])
		assembly = read.csv(file="modules/copy_number/hg19_cytoBandIdeo.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
		plotIdeogram(chrom=ii, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
		close.screen(all.screens=TRUE)
		dev.off()
		
	}
	cat("done!\n", file=paste0("cnvaccess/plot/bychr/", opt$sample_name, "/timestamp"), append=FALSE)

} else if (as.numeric(opt$type)==3) {

	dataA = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-A.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataA[dataA[,"chromosome"]=="X", "chromosome"] = 23
	dataA[dataA[,"chromosome"]=="Y", "chromosome"] = 24
	dataA[,"chromosome"] = as.numeric(dataA[,"chromosome"])
	dataA = subset(dataA, dataA[,"chromosome"]<=23)
	
	tmp = dataA[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataA[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataB = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".pool-B.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataB[dataB[,"chromosome"]=="X", "chromosome"] = 23
	dataB[dataB[,"chromosome"]=="Y", "chromosome"] = 24
	dataB[,"chromosome"] = as.numeric(dataB[,"chromosome"])
	dataB = subset(dataB, dataB[,"chromosome"]<=23)
	
	tmp = dataB[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=2.5, k=10, verbose=FALSE, return.outliers=TRUE)
	dataB[tmp$wins.outliers[,3]!=0,"log2"] = NA
	
	dataC = read.table(file=paste0("cnvaccess/cnr/", opt$sample_name, ".no-pool.cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	dataC[dataC[,"chromosome"]=="X", "chromosome"] = 23
	dataC[dataC[,"chromosome"]=="Y", "chromosome"] = 24
	dataC[,"chromosome"] = as.numeric(dataC[,"chromosome"])
	dataC = subset(dataC, dataC[,"chromosome"]<=23)
	
	tmp = dataC[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=1.5, k=25, verbose=FALSE, return.outliers=TRUE)
	dataC[tmp$wins.outliers[,3]!=0,"log2"] = NA

	CN = rbind(dataA, dataB, dataC)
	index = order(CN[,2])
	CN = CN[index,,drop=FALSE]
	index = order(CN[,1])
	CN = CN[index,,drop=FALSE]
	CN = CN[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
	
	z0 = CN[,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
	tmp = pcf(data=z0, kmin=10, gamma=50, normalize=FALSE, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	save(CN, tmp, file=paste0("cnvaccess/segmented/", opt$sample_name, ".RData"))

} else if (as.numeric(opt$type)==4) {
	
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
	
	'plot_log2_' <- function(x, y, title = "", alpha = NA, psi = NA) {
		par(mar=c(5, 5, 4, 2)+.1)
		data("CytoBand")
		end = NULL
		for (j in 1:23) {
			end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
		}
		end = cumsum(end)
		start = rep(0, 23)
		start[2:23] = end[1:22]+1
		for (j in 1:23) {
			y[y[,"Chromosome"]==j,"Start"] = y[y[,"Chromosome"]==j,"Start"] + start[j]
			y[y[,"Chromosome"]==j,"End"] = y[y[,"Chromosome"]==j,"End"] + start[j]
			x[x[,"Chromosome"]==j,"Position"] = x[x[,"Chromosome"]==j,"Position"] + start[j]
		}
		plot(x[,"Position"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey75", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
		for (j in 1:nrow(y)) {
			lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=1.75, col="red")
		}
		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		abline(v=1, col="goldenrod3", lty=3, lwd=.5)
		abline(h=0, col="red", lty=1, lwd=1)
		for (j in 2:23) {
			v = start[j]
			abline(v=v, col="goldenrod3", lty=3, lwd=.5)
		}
		abline(v=max(x[,"Position"]), col="goldenrod3", lty=3, lwd=.5)
		axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)	
		rect(xleft=1-1e10, xright=x[nrow(x),"Position"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
		title(main = paste0(title, " | alpha = ", signif(alpha, 3), " | psi = ", signif(psi, 3)), line=-1, cex.main=.75, font.main=1)
		box(lwd=1.5)
	}
	
	load(paste0("cnvaccess/segmented/", opt$sample_name, ".RData"))
	tmp = prunesegments.cn(x=tmp, n=5)
	pdf(file=paste0("cnvaccess/plot/segmented/", opt$sample_name, ".pdf"), width=10, height=4.25)
	plot_log2_(x=CN, y=tmp, title = opt$sample_name, alpha=NA, psi=NA)
	dev.off()

}
