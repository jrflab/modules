'plot_sample_lrr_' <- function(x, fit)
{
    mat = x$jointseg
    cncf = fit$cncf
    dipLogR <- fit$dipLogR
    par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
    chr = mat$chrom
    len = table(chr)
    altcol = rep_len(c("light blue", "gray"), length(len))
    chr.col = rep(altcol, len)
    nmark = cncf$num.mark
    tmp = cumsum(len)
    start = c(1, tmp[-length(len)] + 1)
    end = tmp
    mid = start + len/2
    plot(mat$cnlr, pch = ".", axes = F, cex = 1.5, ylim = c(-5,5), col = c("grey", "lightblue")[1 + rep(cncf$chrom - 2 * floor(cncf$chrom/2), cncf$num.mark)],
         ylab = expression(Log[2]~"Ratio"), xlab="Chromosomes")
    points(rep(cncf$cnlr.median, cncf$num.mark), pch = ".", cex = 2, col = "brown")
    labs <- names(mid)
    labs <- sub('21', '', labs)
    labs <- sub('23', 'X', labs)
    axis(side = 1, at = mid, labs, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1, las=2)
    abline(h=0, lty=2, col="lightgrey")
    box()
}

'plot_cncf_' <- function(x, emfit=NULL, clustered=FALSE, plot.type=c("em","naive","both","none"), sname=NULL)
{
    def.par <- par(no.readonly = TRUE)
    plot.type <- match.arg(plot.type)
    if (plot.type=="none") {
    	layout(matrix(1:2, ncol=1))
    }
    if (plot.type=="em") {
    	layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
    }
    if (plot.type=="naive") {
    	layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
    }
    if (plot.type=="both") {
    	layout(matrix(rep(1:6, c(9,9,6,1,6,1)), ncol=1))
    }
    par(mar=c(0.25,3,0.25,1), mgp=c(1.75, 0.6, 0), oma=c(3,0,1.25,0))
    jseg <- x$jointseg
    chrbdry <- which(diff(jseg$chrom) != 0)
    if (missing(emfit)) {
        out <- x$out
        if (plot.type=="em" | plot.type=="both") {
            warning("emfit is missing; plot.type set to naive")
            plot.type <- "naive"
        }
    } else {
        out <- emfit$cncf
        out$tcn <- x$out$tcn
        out$lcn <- x$out$lcn
        out$cf <- x$out$cf
    }
    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    } else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    chrcol <- 1+rep(out$chr-2*floor(out$chr/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0,out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]
    
    plot(jseg$cnlr[is.finite(jseg$cnlr)], pch=".", cex=2, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab=expression(Log[2]~"Ratio"), xaxt="n", las=1, ylim=c(-2,2))
    abline(v=chrbdry, lwd=0.5, lty=3, col="brown")
    abline(v=1, lwd=0.5, lty=3, col="brown")
    abline(v=nrow(jseg$cnlr[is.finite(jseg$cnlr)]), lwd=0.5, lty=3, col="brown")
    abline(h=median(jseg$cnlr, na.rm=TRUE), col="green2")
    abline(h=x$dipLogR, col = "magenta4")
    segments(segstart, cnlr.median, segend, cnlr.median, lwd=1.75, col=2)
  
    plot(jseg$valor[is.finite(jseg$cnlr)], pch=".", cex=2.5, col = c("grey","lightblue","azure4","slateblue")[chrcol], ylab=expression(Log~"OR"), ylim=c(-4,4), xaxt="n", las=1)
    abline(v=chrbdry, lwd=0.5, lty=3, col="brown")
    segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd=1.75, col=2)
    segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd=1.75, col=2)

    # naive copy number and cellular faction pieces
    cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")
    if (plot.type=="naive" | plot.type=="both") {
        # plot the estimated copy numbers and cf
        out$tcn[out$tcn > 10] <- 9 + log10(out$tcn[out$tcn > 10])
        ii <- which(out$lcn > 5)
        if (length(ii)>0) out$lcn[ii] <- 5 + log10(out$lcn[ii])
        plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn)), type="n", ylab="Absolute copies", xaxt="n", las=1)
        abline(v=chrbdry, lwd=0.25)
        segments(segstart, out$lcn, segend, out$lcn, lwd=1.75, col=2)
        segments(segstart, out$tcn, segend, out$tcn, lwd=1.75, col=1)
        # add the cf
        plot(c(0,length(jseg$cnlr)), 0:1, type="n", ylab="CF", xaxt="n", yaxt="n")
        cfcol <- cfpalette[round(10*out$cf+0.501)]
        rect(segstart, 0, segend, 1, col=cfcol, border=NA)
    }
    # EM copy number and cellular faction pieces
    if (plot.type=="em" | plot.type=="both") {
        # plot the estimated copy numbers and cf
        out$tcn.em[out$tcn.em > 10] <- 9 + log10(out$tcn.em[out$tcn.em > 10])
        ii <- which(out$lcn.em > 5)
        if (length(ii)>0) out$lcn.em[ii] <- 5 + log10(out$lcn.em[ii])
        plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn.em)), type="n", ylab="Absolute copies", xaxt="n", las=1)
        abline(v=chrbdry, lwd=0.5, lty=3, col="brown")
        segments(segstart, out$lcn.em, segend, out$lcn.em, lwd=1.75, col=2)
        segments(segstart, out$tcn.em, segend, out$tcn.em, lwd=1.75, col=1)
        # add the cf
        plot(c(0,length(jseg$cnlr)), 0:1, type="n", ylab="CF", xaxt="n", yaxt="n")
        cfcol <- cfpalette[round(10*out$cf.em+0.501)]
        rect(segstart, 0, segend, 1, col=cfcol, border=NA)
    }
    
    # now add the chromosome ticks on x-axis
    chromlevels <- x$chromlevels
    # just make sure chromlevels actually exists
    if (is.null(chromlevels)) chromlevels <- 1:length(nn)
    axis(labels=chromlevels, side=1, at=(nn+c(0,nn[-length(nn)]))/2, cex=0.65)
    # mtext(side=1, line=1.75, "Chromosome", cex=0.8)
    if (!missing(sname)) mtext(sname, side=3, line=0, outer=TRUE, cex=0.8)
    par(def.par)  #- reset to default
}

'plot_log2_' <- function(x, y, n=10, purity=NA, ploidy=NA, title = "")
{

	cn = x$jointseg %>%
		 select(chrom, pos = maploc, log2 = cnlr)
	seg = y$cncf %>%
		  select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	seg = prune_(x=seg, n) %>%
		  mutate(n = cumsum(n))
		  
	purity = ifelse(is.na(purity), 1, purity)
	ploidy = ifelse(is.na(ploidy), 2, ploidy)
	
	data(CytoBand)
   	par(mar=c(5, 5, 4, 2)+.1)
   	end = NULL
   	for (i in 1:23) {
   		end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   	}
   	end = cumsum(end)
   	start = c(1, end[1:22]+1)
   	CytoBand = cbind(start, end)
   	index = NULL
   	for (i in 1:23) {
   		index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(cn$chrom==i)))
   	}
	plot(index, cn$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
 	for (j in 1:nrow(seg)) {
 		if (j == 1) {
 			lines(x=c(1, index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 		} else {
 			lines(x=c(index[seg[j-1,"n"]], index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 		}
  	}
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3", lty=3, lwd=.5)
	for (j in 1:23) {
		abline(v=CytoBand[j,"end"], col="goldenrod3", lty=3, lwd=.5)
	}
	axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	abline(h=0, col="brown", lty=1)
	rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = paste0(title, " | alpha = ", signif(purity, 3), " | psi = ", signif(ploidy, 3)), line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

'psi' <- function (x, z)
{
    xwin <- x
    xwin[x < -z] <- -z
    xwin[x > z] <- z
    return(xwin)
}

'numericChrom' <- function (chrom) 
{
    if (!is.numeric(chrom)) {
        if (is.factor(chrom)) {
            chrom <- as.character(chrom)
        }
        chrx <- c(which(chrom == "x"), which(chrom == "X"))
        chrom[chrx] <- 23
        chry <- c(which(chrom == "y"), which(chrom == "Y"))
        chrom[chry] <- 24
        chrom <- as.numeric(chrom)
    }
    return(chrom)
}

'numericArms' <- function (chrom, char.arms)
{
    p.arm <- which(char.arms == "p")
    q.arm <- which(char.arms == "q")
    arms <- rep(NA, length(char.arms))
    arms[p.arm] <- chrom[p.arm] * 2 - 1
    arms[q.arm] <- chrom[q.arm] * 2
    return(arms)
}

'madWins' <- function (x, tau, k, digits) 
{
    xhat <- medianFilter(x, k)
    d <- x - xhat
    SD <- mad(d)
    z <- tau * SD
    xwin <- xhat + psi(d, z)
    outliers <- rep(0, length(x))
    outliers[round(x, digits) > round(xwin, digits)] <- 1
    outliers[round(x, digits) < round(xwin, digits)] <- -1
    return(list(ywin = xwin, sdev = SD, outliers = outliers))
}

'medianFilter' <- function (x, k) 
{
    n <- length(x)
    filtWidth <- 2 * k + 1
    if (filtWidth > n) {
        if (n == 0) {
            filtWidth <- 1
        }
        else if (n%%2 == 0) {
            filtWidth <- n - 1
        }
        else {
            filtWidth <- n
        }
    }
    runMedian <- runmed(x, k = filtWidth, endrule = "median")
    return(runMedian)
}


'winsorize' <- function (data, pos.unit = "bp", arms = NULL, method = "mad", 
    tau = 2.5, k = 25, gamma = 40, iter = 1, assembly = "hg19", 
    digits = 4, return.outliers = FALSE, save.res = FALSE, file.names = NULL, 
    verbose = TRUE) 
{
    stopifnot(pos.unit %in% c("bp", "kbp", "mbp"))
    stopifnot(method %in% c("mad", "pcf"))
    if (!assembly %in% c("hg19", "hg18", "hg17", "hg16", "mm7", "mm8", "mm9")) {
        stop("assembly must be one of hg19, hg18, hg17 or hg16", call. = FALSE)
    }
    stopifnot(class(data) %in% c("matrix", "data.frame", "character"))
    isfile <- class(data) == "character"
    if (!isfile) {
        stopifnot(ncol(data) >= 3)
        chrom <- data[, 1]
        pos <- data[, 2]
        nSample <- ncol(data) - 2
        sample.names <- colnames(data)[-c(1:2)]
    }
    else {
        f <- file(data, "r")
        head <- scan(f, nlines = 1, what = "character", quiet = TRUE, 
            sep = "\t")
        if (length(head) < 3) {
            stop("Data in file must have at least 3 columns", 
                call. = FALSE)
        }
        sample.names <- head[-c(1:2)]
        nSample <- length(sample.names)
        chrom.pos <- read.table(file = data, sep = "\t", header = TRUE, 
            colClasses = c(rep(NA, 2), rep("NULL", nSample)), 
            as.is = TRUE)
        chrom <- chrom.pos[, 1]
        pos <- chrom.pos[, 2]
    }
    if (is.factor(chrom)) {
        chrom <- as.character(chrom)
    }
    num.chrom <- numericChrom(chrom)
    nProbe <- length(num.chrom)
    if (!is.numeric(pos)) {
        stop("input in data column 2 (posistions) must be numeric", call. = FALSE)
    }
    if (is.null(arms)) {
        arms <- getArms(num.chrom, pos, pos.unit, get(assembly))
    }
    else {
        if (length(arms) != nProbe) {
            stop("'arms' must be the same length as number of rows in data", 
                call. = FALSE)
        }
    }
    num.arms <- numericArms(num.chrom, arms)
    arm.list <- unique(num.arms)
    nArm <- length(arm.list)
    if (!save.res) {
        wins.data <- matrix(nrow = 0, ncol = nSample)
        if (return.outliers) {
            wins.outliers <- matrix(nrow = 0, ncol = nSample)
        }
    }
    else {
        if (is.null(file.names)) {
            dir.res <- "Wins_res"
            if (!dir.res %in% dir()) {
                dir.create(dir.res)
            }
            file.names <- c(paste(dir.res, "/", "wins.data.txt", 
                sep = ""), paste(dir.res, "/", "wins.outliers.txt", 
                sep = ""))
        }
        else {
            if (length(file.names) < 2) {
                stop("'file.names' must be of length 2", call. = FALSE)
            }
        }
    }
    for (c in 1:nArm) {
        probe.c <- which(num.arms == arm.list[c])
        wins.data.c <- matrix(nrow = length(probe.c), ncol = 0)
        if (return.outliers || save.res) {
            wins.outliers.c <- matrix(nrow = length(probe.c), 
                ncol = 0)
        }
        if (!isfile) {
            arm.data <- data[probe.c, -c(1:2), drop = FALSE]
        }
        else {
            arm.data <- read.table(f, nrows = length(probe.c), 
                sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", 
                  nSample)))
        }
        if (any(!sapply(arm.data, is.numeric))) {
            stop("input in data columns 3 and onwards (copy numbers) must be numeric", 
                call. = FALSE)
        }
        for (i in 1:nSample) {
            y <- arm.data[, i]
            na <- is.na(y)
            use.y <- y[!na]
            ywins <- rep(NA, length(y))
            outliers <- rep(NA, length(y))
            wins <- switch(method, mad = madWins(use.y, tau = tau, 
                k = k, digits = digits), pcf = pcfWins(use.y, 
                tau = tau, k = k, gamma = gamma, iter = iter, 
                digits = digits))
            ywins[!na] <- wins$ywin
            outliers[!na] <- wins$outliers
            ywins <- round(ywins, digits = digits)
            wins.data.c <- cbind(wins.data.c, ywins)
            if (return.outliers || save.res) {
                wins.outliers.c <- cbind(wins.outliers.c, outliers)
            }
        }
        if (!save.res) {
            wins.data <- rbind(wins.data, wins.data.c)
            if (return.outliers) {
                wins.outliers <- rbind(wins.outliers, wins.outliers.c)
            }
        }
        else {
            if (c == 1) {
                wd <- file(file.names[1], "w")
                wo <- file(file.names[2], "w")
            }
            write.table(data.frame(chrom[probe.c], pos[probe.c], 
                wins.data.c, stringsAsFactors = FALSE), file = wd, 
                col.names = if (c == 1) 
                  c("chrom", "pos", sample.names)
                else FALSE, row.names = FALSE, quote = FALSE, 
                sep = "\t")
            write.table(data.frame(chrom[probe.c], pos[probe.c], 
                wins.outliers.c, stringsAsFactors = FALSE), file = wo, 
                col.names = if (c == 1) 
                  c("chrom", "pos", sample.names)
                else FALSE, row.names = FALSE, quote = FALSE, 
                sep = "\t")
        }
        if (verbose) {
            chr <- unique(chrom[probe.c])
            a <- unique(arms[probe.c])
            cat(paste("winsorize finished for chromosome arm ", 
                chr, a, sep = ""), "\n")
        }
    }
    if (isfile) {
        close(f)
    }
    if (!save.res) {
        wins.data <- data.frame(chrom, pos, wins.data, stringsAsFactors = FALSE)
        colnames(wins.data) <- c("chrom", "pos", sample.names)
        if (return.outliers) {
            wins.outliers <- data.frame(chrom, pos, wins.outliers, 
                stringsAsFactors = FALSE)
            colnames(wins.outliers) <- c("chrom", "pos", sample.names)
            return(list(wins.data = wins.data, wins.outliers = wins.outliers))
        }
        else {
            return(wins.data)
        }
    }
    else {
        close(wd)
        close(wo)
        cat(paste("winsorized data were saved in file", file.names[1]), 
            sep = "\n")
        cat(paste("outliers were saved in file", file.names[2]), 
            sep = "\n")
        return(invisible(NULL))
    }
}

'prune_' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"log2"] - 2^x[,"log2"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
 			mcl = mean(x[indx,"log2"])
			scl = sd(x[indx,"log2"])
			ind = which(x[indx,"log2"]<(mcl+1.96*scl) & x[indx,"log2"]>(mcl-1.96*scl))
			x[indx[ind],"log2"] = mean(x[indx[ind],"log2"])
		} else {
			x[indx,"log2"] = mean(x[indx,"log2"])
		}
	}
	return(x)
}
