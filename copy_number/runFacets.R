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

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
                make_option("--pre_cval", default = 50, type = 'integer', help = "pre-processing critical value"),
                make_option("--cval", default = 50, type = 'integer', help = "critical value for segmentation"),
                make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
                make_option("--genome", default = 'hg19', type = 'character', help = "genome of counts file"),
                make_option("--outPrefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

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


out <- baseCountFile %>% preProcSample(snp.nbhd = opt$snp_nbhd, cval = opt$pre_cval, chromlevels = chromLevels) %>%
    procSample(cval = opt$cval, min.nhet = opt$min_nhet)
fit <- out %>% emcncf

CairoPNG(file = str_c(opt$outPrefix,".biseg.png"), height = 1000, width = 800)
plotSample(out, chromlevels = chromLevels)
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
out$IGV = formatSegmentOutput(out, id)
save(out, fit, file = str_c(opt$outPrefix, ".Rdata"), compress=T)

ff = str_c(opt$outPrefix, ".out")
cat("# Version =", version, "\n", file = ff, append = T)
cat("# Input =", basename(baseCountFile), "\n", file = ff, append = T)
cat("# tumor =", tumorName, "\n", file = ff, append = T)
cat("# normal =", normalName, "\n", file = ff, append = T)
cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = T)
cat("# cval =", opt$cval, "\n", file = ff, append = T)
cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = T)
cat("# genome =", opt$genome, "\n", file = ff, append = T)
cat("# Purity =", fit$purity, "\n", file = ff, append = T)
cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = T)
cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = T)
cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
cat("# loglik =", fit$loglik, "\n", file = ff, append = T)

write.table(cbind(out$IGV[, 1:4], fit$cncf[, 2:ncol(fit$cncf)]), 
    str_c(opt$outPrefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')

CairoPNG(file = str_c(opt$outPrefix, ".cncf.png"), height = 1100, width = 850)
plotSampleCNCF(out, fit)
#plotSampleCNCF.custom(out$jointseg, out$out, fit, 
#        main = paste(projectName, "[", tumorName, normalName, "]", "cval  = ", CVAL))

dev.off()
warnings()

