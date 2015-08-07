#!/usr/bin/env Rscript
# run the facets library

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
suppressPackageStartupMessages(library("Cairo"));

plotSampleCNCF <- function (x, fit) 
{
    mat = x$jointseg
    cncf = fit$cncf
    dipLogR <- fit$dipLogR
    layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), ncol = 1))
    par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
    chr = mat$chrom
    len = table(chr)
    altcol = rep(c("light blue", "gray"), 12)[-24]
    chr.col = rep(altcol, len)
    nmark = cncf$num.mark
    tmp = cumsum(len)
    start = c(1, tmp[-23] + 1)
    end = tmp
    mid = start + len/2
    plot(mat$cnlr, pch = ".", axes = F, cex = 1.5, ylim = c(-3, 
        3), col = c("grey", "lightblue")[1 + rep(cncf$chrom - 
        2 * floor(cncf$chrom/2), cncf$num.mark)], ylab = "log-ratio")
    abline(h = dipLogR, col = "magenta4")
    points(rep(cncf$cnlr.median, cncf$num.mark), pch = ".", cex = 2, 
        col = "brown")
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1)
    box()
    plot(mat$valor, axes = F, pch = ".", cex = 1.5, col = c("grey", 
        "lightblue")[1 + rep(cncf$chrom - 2 * floor(cncf$chrom/2), 
        cncf$num.mark)], ylab = "log-odds-ratio", ylim = c(-4, 
        4))
    points(rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch = ".", 
        cex = 2, col = "brown")
    points(-rep(sqrt(abs(cncf$mafR)), cncf$num.mark), pch = ".", 
        cex = 2, col = "brown")
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1)
    box()
    plot(rep(cncf$cf.em, cncf$num.mark), axes = F, pch = ".", 
        cex = 2, xlab = "Chromosome", ylab = "Cellular fraction (EM)", 
        ylim = c(0, 1))
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    tcnscaled <- cncf$tcn.em
    tcnscaled[cncf$tcn.em > 5 & !is.na(cncf$tcn.em)] = (5 + (tcnscaled[cncf$tcn.em > 
        5 & !is.na(cncf$tcn.em)] - 5)/3)
    matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn.em, 
        cncf$num.mark) - 0.1), pch = ".", cex = 3, col = 1:2, 
        lwd = 1, ylab = "Integer copy number (EM)", yaxt = "n", 
        xaxt = "n")
    axis(2, at = c(0:5, 5 + (1:35)/3), labels = 0:40, cex.axis = 1)
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    abline(h = c(0:5, 5 + (1:35)/3), lty = 3, col = "gray")
    plot(rep(cncf$cf, cncf$num.mark), axes = F, pch = ".", cex = 2, 
        xlab = "Chromosome", ylab = "Cellular fraction (cncf)", 
        ylim = c(0, 1))
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    axis(side = 2, cex.axis = 1)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    tcnscaled <- cncf$tcn
    tcnscaled[cncf$tcn > 5 & !is.na(cncf$tcn)] = (5 + (tcnscaled[cncf$tcn > 
        5 & !is.na(cncf$tcn)] - 5)/3)
    matplot(cbind(rep(tcnscaled, cncf$num.mark), rep(cncf$lcn, 
        cncf$num.mark) - 0.1), pch = ".", cex = 3, col = 1:2, 
        lwd = 1, ylab = "Integer copy number (cncf)", yaxt = "n", 
        xaxt = "n")
    axis(2, at = c(0:5, 5 + (1:35)/3), labels = 0:40, cex.axis = 1)
    axis(side = 1, at = mid, 1:23, cex.axis = 1, las = 2)
    box()
    abline(v = start, lty = 3, col = "gray")
    abline(v = end, lty = 3, col = "gray")
    abline(h = c(0:5, 5 + (1:35)/3), lty = 3, col = "gray")
}
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
                make_option("--pre_cval", default = 50, type = 'integer', help = "pre-processing critical value"),
                make_option("--cval1", default = 150, type = 'integer', help = "critical value for estimating diploid log Ratio"),
                make_option("--cval2", default = 50, type = 'integer', help = "starting critical value for segmentation (increases by 10 until success)"),
                make_option("--maxCval", default = 300, type = 'integer', help = "maximum critical value for segmentation (increases by 10 until success)"),
                make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
                make_option("--genome", default = 'hg19', type = 'character', help = "genome of counts file"),
                make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);
n

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need base counts file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    baseCountFile <- arguments$args[1];
}

tumorName <- baseCountFile %>% sub('.*/', '', .) %>% sub('_.*', '', .)
normalName <- baseCountFile %>% sub('.*/', '', .) %>% sub('.*_', '', .) %>% sub('\\..*', '', .)

switch(opt$genome,
       hg19={
           data(hg19gcpct)
           chromLevels=c(1:22, "X")
       },
       mm9={
           data(mm9gcpct)
           chromLevels=c(1:19)
       },
       {
           stop(paste("Invalid Genome",opt$genome))
       })



buildData=installed.packages()["facets",]
cat("#Module Info\n")
for(fi in c("Package","LibPath","Version","Built")){
    cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
}
version=buildData["Version"]
cat("\n")


preOut <- baseCountFile %>% preProcSample(snp.nbhd = opt$snp_nbhd, cval = opt$pre_cval, chromlevels = chromLevels)
out1 <- preOut %>% procSample(cval = opt$cval1, min.nhet = opt$min_nhet)

cval <- opt$cval2
success <- F
while (!success || cval > opt$maxCval) {
    out2 <- preOut %>% procSample(cval = cval, min.nhet = opt$min_nhet, dipLogR = out1$dipLogR)
    print(str_c("attempting to run emncf() with cval = ", cval))
    fit <- tryCatch({
        out2 %>% emcncf
    }, error = function(e) {
        print(paste("Error:", e))
        return(NULL)
    })
    if (!is.null(fit)) {
        success <- T
    } else {
        cval <- cval + 10
    }
}
if (!success) {
    stop("Failed to segment data\n")
}


CairoPNG(file = str_c(opt$outPrefix,".biseg.png"), height = 1000, width = 800)
plotSample(out2, chromlevels = chromLevels)
dev.off()

pdf(file = str_c(opt$outPrefix, ".biseg.pdf"), height = 12, width = 9)
plotSample(out2, chromlevels = chromLevels)
dev.off()

formatSegmentOutput=function(out,sampID) {
	seg=list()
	seg$ID=rep(sampID,nrow(out$out))
	seg$chrom=out$out$chr
	seg$loc.start=rep(NA,length(seg$ID))
	seg$loc.end=seg$loc.start
	seg$num.mark=out$out$num.mark
	seg$seg.mean=out$out$cnlr.median
	for(i in 1:nrow(out$out)) {
		lims=range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chr[i] & out$jointseg$seg==out$out$seg[i])],na.rm=T)
		seg$loc.start[i]=lims[1]
		seg$loc.end[i]=lims[2]
	}
	as.data.frame(seg)
}
id <- paste(tumorName, normalName, sep = '_')
out2$IGV = formatSegmentOutput(out2, id)
save(out2, fit, file = str_c(opt$outPrefix, ".Rdata"), compress=T)

ff = str_c(opt$outPrefix, ".out")
cat("# Version =", version, "\n", file = ff, append = T)
cat("# Input =", basename(baseCountFile), "\n", file = ff, append = T)
cat("# tumor =", tumorName, "\n", file = ff, append = T)
cat("# normal =", normalName, "\n", file = ff, append = T)
cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = T)
cat("# cval1 =", opt$cval1, "\n", file = ff, append = T)
cat("# cval2 =", cval, "\n", file = ff, append = T)
cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = T)
cat("# genome =", opt$genome, "\n", file = ff, append = T)
cat("# Purity =", fit$purity, "\n", file = ff, append = T)
cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = T)
cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = T)
cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
cat("# loglik =", fit$loglik, "\n", file = ff, append = T)

write.table(cbind(out2$IGV[, 1:4], fit$cncf[, 2:ncol(fit$cncf)]), 
    str_c(opt$outPrefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')

CairoPNG(file = str_c(opt$outPrefix, ".cncf.png"), height = 1100, width = 850)
plotSampleCNCF(out2, fit)
dev.off()

pdf(file = str_c(opt$outPrefix, ".cncf.pdf"), height = 12, width = 7)
plotSampleCNCF(out2, fit)
dev.off()

#plotSampleCNCF.custom(out$jointseg, out$out, fit, 
#        main = paste(projectName, "[", tumorName, normalName, "]", "cval  = ", CVAL))
warnings()

