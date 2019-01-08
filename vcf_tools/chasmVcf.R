#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
}
#options(error = recover)
#options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'b37', help = "genome build [default %default]"),
                make_option("--chasmDir", default = NULL, help = "CHASM dir"),
                make_option("--snvBoxDir", default = NULL, help = "snvBox dir"),
                make_option("--classifier", default = 'Breast', help = "CHASM classifier(s), comma-delimited [default %default]"),
                make_option("--outFile", default = NULL, help = "vcf output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$chasmDir)) {
    cat("Need CHASM dir\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

classifiers <- unlist(strsplit(opt$classifier, ','))

# create new header
vcfHeader <- scanVcfHeader(fn)
hinfoprime <- as.data.frame(info(vcfHeader))
for (cl in classifiers) {
    cln <- gsub('-', '_', cl)
    X <- rbind(chasm_mut = c("Number" = "A", "Type" = "String", "Description" = paste(cln, "CHASM mutation")),
               chasm_pval = c("Number" = "A", "Type" = "Float", "Description" = paste(cln, "CHASM p-value")),
               chasm_score = c("Number" = "A", "Type" = "Float", "Description" = paste(cln, "CHASM score")),
               chasm_pred = c("Number" = "A", "Type" = "String", "Description" = paste(cln, "CHASM prediction")),
               chasm_fdr = c("Number" = "A", "Type" = "Float", "Description" = paste(cln, "CHASM B-H FDR")))
    rownames(X) <- paste(cln, rownames(X), sep = "_")
    hinfoprime <- rbind(hinfoprime, X)
}
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(vcfHeader)
hlist$INFO <- hinfoprime
colnames(hlist$INFO) <- c('Number', 'Type', 'Description')
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)


temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
tab <- TabixFile(zipped, idx, yieldSize = 8000)
open(tab)
while(nrow(vcf <- readVcf(tab, genome = opt$genome))) {
    metadata(vcf)$header <- newVcfHeader

    oldwd <- getwd()

    if (sum(rowRanges(vcf)$FILTER == "PASS") > 0) {
        # convert vcf to chasm input
        al <- sapply(rowRanges(vcf)$ALT, length)
        df <- data.frame(id = 1:nrow(vcf), seq = as.character(seqnames(rowRanges(vcf))), zero = as.integer(start(rowRanges(vcf)) - 1), one = as.integer(start(rowRanges(vcf))), strand = rep('+', nrow(vcf)), ref = as.character(rowRanges(vcf)$REF), filt = rowRanges(vcf)$FILTER, stringsAsFactors = F)
        re <- rep.int(df$id, times = al)
        df <- df[re, ,drop = F]
        df$alt <- as.character(unlist(rowRanges(vcf)$ALT))
        df <- subset(df, sapply(alt, nchar) == 1 & sapply(ref, nchar) == 1 & filt == "PASS")
        df <- df[,-which(colnames(df) == 'filt')]
        if (all(!grepl('chr', df$seq))) {
            df$seq <- sub('^', 'chr', df$seq)
        }

        setwd(opt$chasmDir)
        tmp <- tempfile()
        co <- c('id', 'seq', 'one', 'strand', 'ref', 'alt')
        write.table(df[, co], file = tmp, quote = F, sep = '\t', row.names = F, col.names = F)
        for (cl in classifiers) {
            cln <- gsub('-', '_', cl)
            cmd <- paste("CHASMDIR=", opt$chasmDir, ' SNVBOXDIR=', opt$snvBoxDir, ' python RunChasm ', cl, ' ', tmp, ' -g' ,sep = '')
            #cat(cmd)
            system(cmd, ignore.stdout = T)
            results <- read.table(file = paste(tmp, cl, '.output', sep = ''), sep = '\t', header = T, as.is = T)
            if (nrow(results) > 1) {
                infoprime <- info(vcf)
                infoprime[as.integer(as.character(results$MutationID)), paste(cln, "chasm_mut", sep = "_")] <- as.character(results$Mutation)
                infoprime[as.integer(as.character(results$MutationID)), paste(cln, "chasm_score", sep = "_")] <- results$CHASM
                infoprime[as.integer(as.character(results$MutationID)), paste(cln, "chasm_pred", sep = "_")] <- ifelse(results$CHASM <= 0.3, 'Driver', 'Passenger')
                infoprime[as.integer(as.character(results$MutationID)), paste(cln, "chasm_pval", sep = "_")] <- results$PValue
                if (nrow(results) > 10) {
                    infoprime[as.integer(as.character(results$MutationID)), paste(cln, "chasm_fdr", sep = "_")] <- results$BHFDR
                }
                info(vcf) <- infoprime
            }
        }
    }

    # fix sample genotype order
    if ("GT" %in% names(geno(vcf))) {
        x <- which(names(geno(vcf)) == "GT")
        ord <- c(x, (1:length(geno(vcf)))[-x])
        geno(vcf) <- geno(vcf)[ord]
    }

    setwd(oldwd)
    writeVcf(vcf, out)
}
close(tab)
close(out)
