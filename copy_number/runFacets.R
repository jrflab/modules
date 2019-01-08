#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("facets"))
suppressPackageStartupMessages(library("pctGCdata"))
suppressPackageStartupMessages(library("foreach"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
                make_option("--pre_cval", default = 50, type = 'integer', help = "pre-processing critical value"),
                make_option("--cval1", default = 150, type = 'integer', help = "critical value for estimating diploid log Ratio"),
                make_option("--cval2", default = 50, type = 'integer', help = "starting critical value for segmentation (increases by 10 until success)"),
                make_option("--max_cval", default = 5000, type = 'integer', help = "maximum critical value for segmentation (increases by 10 until success)"),
                make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
                make_option("--het_threshold", default = 0.25, type = 'double', help = "AF threshold for heterozygous SNPs"),
                make_option("--diplogr", default = NULL, type = 'double', help = "override diploid log-ratio"),
                make_option("--ndepth_max", default = 1000, type = 'integer', help = "normal depth max"),
                make_option("--use_emcncf2", default = F, action = 'store_true', help = "use emcncf version 2"),
                make_option("--gene_loc_file", default = '~/share/reference/IMPACT410_genes_for_copynumber.txt', type = 'character', help = "file containing gene locations"),
                make_option("--genome", default = 'b37', type = 'character', help = "genome of counts file"),
                make_option("--out_prefix", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need base counts file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$out_prefix)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    snpPileupFile <- arguments$args[1];
}

tumorName <- snpPileupFile %>% sub('.*/', '', .) %>% sub('_.*', '', .)
normalName <- snpPileupFile %>% sub('.*/', '', .) %>% sub('^.*_', '', .) %>% sub('\\..*', '', .)

switch(opt$genome,
       b37={
            facetsGenome = 'hg19'
       },
       GRCh37={
            facetsGenome = 'hg19'
       },
       hg19={
            facetsGenome = 'hg19'
       },
       mm9={
            facetsGenome = 'mm9'
       },
       mm10={
            facetsGenome = 'mm10'
       },
       GRCm38={
            facetsGenome = 'mm10'
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


snpmat <- readSnpMatrix(snpPileupFile)
preOut <- snpmat %>% preProcSample(snp.nbhd = opt$snp_nbhd, het.thresh = opt$het_threshold, cval = opt$pre_cval, gbuild = facetsGenome, ndepthmax = opt$ndepth_max)
if (!is.null(opt$diplogr)) {
    cval <- opt$cval2
    out2 <- preOut %>% procSample(cval = cval, min.nhet = opt$min_nhet, dipLogR = opt$diplogr)
    if (opt$use_emcncf2) {
        fit <- out2 %>% emcncf2
    } else {
        fit <- out2 %>% emcncf
    }
} else {
    out1 <- preOut %>% procSample(cval = opt$cval1, min.nhet = opt$min_nhet)

    cval <- opt$cval2
    success <- F
    while (!success && cval < opt$max_cval) {
        out2 <- preOut %>% procSample(cval = cval, min.nhet = opt$min_nhet, dipLogR = out1$dipLogR)
        print(str_c("attempting to run emncf() with cval2 = ", cval))
        fit <- tryCatch({
            if (opt$use_emcncf2) {
                out2 %>% emcncf2
            } else {
                out2 %>% emcncf
            }
        }, error = function(e) {
            print(paste("Error:", e))
            return(NULL)
        })
        if (!is.null(fit)) {
            success <- T
        } else {
            cval <- cval + 100
        }
    }
    if (!success) {
        stop("Failed to segment data\n")
    }
}

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
save(out2, fit, file = str_c(opt$out_prefix, ".Rdata"), compress=T)

ff = str_c(opt$out_prefix, ".out")
cat("# Version =", version, "\n", file = ff, append = T)
cat("# Input =", basename(snpPileupFile), "\n", file = ff, append = T)
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

tab <- cbind(select(out2$IGV, ID:num.mark), select(fit$cncf, -start, -end, -chrom, -num.mark))
write.table(tab, str_c(opt$out_prefix, ".cncf.txt"), row.names = F, quote = F, sep = '\t')

write.table(out2$IGV, str_c(opt$out_prefix, '.facets.seg'), row.names = F, quote = F, sep = '\t')

warnings()

