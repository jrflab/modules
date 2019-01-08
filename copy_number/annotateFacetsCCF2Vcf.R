#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("VariantAnnotation"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--ccfRscript", default = '~/share/usr/ccf.R', help='computCCF and confCCF R script'),
                make_option("--genome", default = 'b37', type = 'character', help = "genome of counts file"),
                make_option("--tumor", default = 'TUMOR', type = 'character', help = "tumor sample"),
                make_option("--purity", default = NULL, type = 'float', help = "purity of sample if overriding facets purity"),
                make_option("--facetsRdata", default = NULL, type = "character", action = "store", help ="facets Rdata file"),
                make_option("--outFile", default = NULL, help = "output prefix"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need vcf file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$facetsRdata)) {
    cat("Need facets Rdata file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n\n")
    print_help(parser);
    stop();
}

source(opt$ccfRscript)

tumorSample <- opt$tumor
vcfFile <- arguments$args[1];
vcf <- readVcf(vcfFile, opt$genome)

# read primary rdata file
load(opt$facetsRdata)
facetsSeg <- fit$cncf
facetsSeg$chrom <- as.character(facetsSeg$chrom)
facetsSeg$chrom[facetsSeg$chrom == '23'] <- 'X'
if (any(grepl('chr', as.character(seqnames(vcf))))) {
    chr <- paste('chr', facetsSeg$chrom, sep = '')
} else {
    chr <- facetsSeg$chrom
}
facetsGr <- with(facetsSeg, GRanges(seqnames = chr,
                                    ranges = IRanges(start = start, end = end),
                                    cf = cf.em, tcn = tcn.em, lcn = lcn.em, mafR = mafR))
medianMafR <- median(facetsGr$mafR)
sdMafR <- sd(facetsGr$mafR)
if (!is.null(opt$purity)) {
    purity <- opt$purity
} else {
    purity <- fit$purity
}


info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "facets cellular fraction",
                                                        row.names = "facetsCF"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "facets log ratio",
                                                        row.names = "facetsMafR"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets total copy number",
                                                        row.names = "facetsTCN_EM"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets lesser copy number",
                                                        row.names = "facetsLCN_EM"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String",
                                                        Description = "facets call",
                                                        row.names = "facetsLOHCall"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "0", Type = "Flag",
                                                        Description = "facets LOH",
                                                        row.names = "facetsLOH"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "String",
                                                        Description = "clonal status",
                                                        row.names = "clonalStatus"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "CCF confidence interval upper bound",
                                                        row.names = "ccfConfUpper"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "CCF confidence interval lower bound",
                                                        row.names = "ccfConfLower"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Float",
                                                        Description = "ccf", row.names = "ccf"))
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number = "1", Type = "Integer",
                                                        Description = "facets multiplicity", row.names = "facetsMultiplicity"))

ref <- sapply(geno(vcf)$AD[, tumorSample], function(x) x[1])
alt <- sapply(geno(vcf)$AD[, tumorSample], function(x) x[2])
vaf <- alt / (alt + ref)
pass <- rowRanges(vcf)$FILTER == 'PASS' & !is.na(vaf)
if (sum(pass) == 0) {
    cat("No unfiltered variants\n")
} else {
    cat(sum(pass), "variants pass\n")
    ol <- findOverlaps(rowRanges(vcf[pass]), facetsGr, select = 'first')
    ol[is.na(ol)] <- nearest(rowRanges(vcf[pass])[is.na(ol)], facetsGr) # get nearest segment for missing values
    if (sum(!is.na(ol)) > 0) {
        tcn <- facetsGr[ol[!is.na(ol)], ]$tcn
        lcn <- facetsGr[ol[!is.na(ol)], ]$lcn
        purity <- rep(purity, length(tcn))

        ref <- sapply(geno(vcf[pass])$AD[!is.na(ol), tumorSample], function(x) x[1])
        alt <- sapply(geno(vcf[pass])$AD[!is.na(ol), tumorSample], function(x) x[2])
        vaf <- alt / (alt + ref)

        ccfFit <- computeCCF(vaf = vaf, tcn, lcn, purity = purity)
        conf <- confCCF(alt = alt, ref = ref, tcn, lcn, purity = purity,
                            multiplicity = ccfFit$multiplicity)
        ccfLower <- conf$lower
        ccfUpper <- conf$upper
        clonalStatus <- ifelse(round(ccfLower, 2) >= 0.75, "clonal", 
                            ifelse(round(ccfLower, 2) < 0.75 & ccfFit$ccf >= 0.8, 'likely_clonal', 
                                    "subclonal"))

        info(vcf)$facetsCF[pass][!is.na(ol)] <- facetsGr$cf[ol[!is.na(ol)]]
        info(vcf)$facetsTCN_EM[pass][!is.na(ol)] <- facetsGr$tcn[ol[!is.na(ol)]]
        info(vcf)$facetsLCN_EM[pass][!is.na(ol)] <- facetsGr$lcn[ol[!is.na(ol)]]
        info(vcf)$facetsMafR[pass][!is.na(ol)] <- facetsGr$mafR[ol[!is.na(ol)]]
        info(vcf)$facetsLOH[pass][!is.na(ol)] <- facetsGr$lcn[ol[!is.na(ol)]] == 0
        info(vcf)$facetsLOHCall[pass][!is.na(ol)] <- ifelse(facetsGr$lcn[ol[!is.na(ol)]] == 0, 'true', 'false')
        info(vcf)$facetsMultiplicity[pass][!is.na(ol)] <- ccfFit$multiplicity
        info(vcf)$ccf[pass][!is.na(ol)] <- ccfFit$ccf
        info(vcf)$clonalStatus[pass][!is.na(ol)] <- clonalStatus
        info(vcf)$ccfConfUpper[pass][!is.na(ol)] <- ccfUpper
        info(vcf)$ccfConfLower[pass][!is.na(ol)] <- ccfLower
        x <- is.na(info(vcf)$facetsLOH[pass][!is.na(ol)])
        if (sum(x) > 0) {
            info(vcf)$facetsLOH[pass][!is.na(ol)][x] <- facetsGr$mafR[ol[!is.na(ol)]][x] > medianMafR + sdMafR
            info(vcf)$facetsLOHCall[pass][!is.na(ol)][x] <- ifelse(facetsGr$mafR[ol[!is.na(ol)]][x] > medianMafR + sdMafR, 'true', 'false')
        }
    }
}


writeVcf(vcf, opt$outFile)
