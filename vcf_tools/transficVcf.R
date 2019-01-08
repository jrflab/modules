#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
        make_option("--transfic", default = "transf_scores.pl", help = "transfic perl script"),
        make_option(c('-g', "--grouping"), default = "cp", help = 'transfic grouping of genes used to compute baseline tolerance [default %default] [options gosmf gosbp cp doms]'),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]"))
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}


fn <- arguments$args[1];
cat("Reading vcf file", fn, "...")
vcf <- readVcf(fn, opt$genome)
cat("done\n")

cols <- c("dbNSFP_Ensembl_geneid", "dbNSFP_SIFT_score", "dbNSFP_Polyphen2_HVAR_score", "dbNSFP_MutationAssessor_score")
passIds <- which(apply(as.data.frame(info(vcf)[,cols]), 1, function(x) all(unlist(x) != ".") && all(!is.na(unlist(x)))))

results <- NULL
if (length(passIds) > 0) {
    tmp1 <- tempfile()
    tmp2 <- tempfile()
    X <- as.data.frame(info(vcf)[passIds, cols])
    X[,1] <- sapply(X[,1], function(x) x[1])
    X[,1] <- sapply(X[,1], function(x) unlist(strsplit(x, "\\|"))[1])
    X[,2] <- sapply(X[,2], function(x) min(as.numeric(unlist(x))))
    if (any(grepl("\\|", X[,3]))) {
        X[,3] <- sapply(X[,3], function(x) max(as.numeric(unlist(strsplit(x, "\\|")))))
    } else {
        X[,3] <- sapply(X[,3], function(x) max(as.numeric(unlist(x))))
    }
    X[,4] <- sapply(X[,4], function(x) max(as.numeric(unlist(x))))
    write.table(X, file = tmp1, quote = F, sep = '\t', col.names = F, row.names = T)
    cmd <- paste(opt$transfic, opt$grouping, tmp1, ">", tmp2)
    system(cmd)
    results <- read.table(tmp2, sep = '\t', as.is = T)
    newCols <- c('id', 'geneId', 'transficSIFT_score', 'transficSIFT_pred', 'transficPolyphen2_HVAR_score', 'transficPolyphen2_HVAR_pred', 'transficMutationAssessor_score', 'transficMutationAssessor_pred')
    mergeCols <- newCols[3:length(newCols)]
    colnames(results) <- newCols
}



hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
rownames(hinfoprime) <- rownames(info(header(vcf)))
hinfoprime <- rbind(hinfoprime, transficSIFT_score = c("A", "Float", "transfic SIFT score"))
hinfoprime <- rbind(hinfoprime, transficPolyphen2_HVAR_score = c("A", "Float", "transfic polyphen2 SIFT score"))
hinfoprime <- rbind(hinfoprime, transficMutationAssessor_score = c("A", "Float", "MutationAssessor score"))
hinfoprime <- rbind(hinfoprime, transficSIFT_pred = c("A", "String", "transfic SIFT prediction"))
hinfoprime <- rbind(hinfoprime, transficPolyphen2_HVAR_pred = c("A", "String", "transfic polyphen2 SIFT prediction"))
hinfoprime <- rbind(hinfoprime, transficMutationAssessor_pred = c("A", "String", "MutationAssessor prediction"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(metadata(vcf)$header)
hlist$INFO <- hinfoprime
metadata(vcf)$header <- new("VCFHeader", samples = header(vcf)@samples, header = hlist)

newScores <- c('transficSIFT_score', 'transficPolyphen2_HVAR_score', 'transficMutationAssessor_score')
newPreds <- c('transficSIFT_pred', 'transficPolyphen2_HVAR_pred', 'transficMutationAssessor_pred')

infoprime <- info(vcf)
for (sc in newScores) {
    infoprime[, sc] <- as.numeric(rep(NA, nrow(vcf)))
}
for (pred in newPreds) {    
    infoprime[, pred] <- as.character(rep(NA, nrow(vcf)))
}
info(vcf) <- infoprime

if (!is.null(results) && nrow(results) > 0) {
    cat("Merging transfic results ... ")
    infodprime <- info(vcf[passIds, ])[results[,1], ]
    for (co in mergeCols) {
        infodprime[, co] <- results[, co]
    }
    info(vcf)[passIds[results[,1]], ] <- infodprime
    cat("done\n")
} else {
    cat("No results from transfic\n")
}

#fix sample genotype order

if ("GT" %in% names(geno(vcf))) {
    x <- which(names(geno(vcf)) == "GT")
    ord <- c(x, (1:length(geno(vcf)))[-x])
    geno(vcf) <- geno(vcf)[ord]
}

cat("Writing vcf to", opt$outFile, "... ")
writeVcf(vcf, opt$outFile)
cat("done\n")


