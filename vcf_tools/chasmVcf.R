#!/usr/bin/env Rscript
# Read a variant table and extract uniprot accession ids

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
#options(error = recover)
#options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--chasmDir", default = NULL, help = "CHASM dir"),
                make_option("--python", default = NULL, help = "python binary"),
                make_option("--classifier", default = 'Breast', help = "CHASM classifier(s), comma-delimited [default %default]"),
                make_option("--outFile", default = NULL, help = "vcf output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$chasmDir)) {
    cat("Need CHASM dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$python)) {
    cat("Need python\n");
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

classifiers <- unlist(strsplit(',', opt$classifier))

# create new header
vcfHeader <- scanVcfHeader(fn)
hinfoprime <- apply(as.data.frame(info(vcfHeader)), 2, as.character)
rownames(hinfoprime) <- rownames(info(vcfHeader))
for (cl in classifiers) {
    X <- rbind(chasm_mut = c("A", "String", "CHASM mutation"),
               chasm_pval = c("A", "Float", "CHASM p-value"),
               chasm_score = c("A", "Float", "CHASM score"),
               chasm_fdr = c("A", "Float", "CHASM B-H FDR"))
    rownames(X) <- paste(cl, rownames(X), sep = "_")
    hinfoprime <- rbind(hinfoprime, X)
}
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(vcfHeader)
hlist$INFO <- hinfoprime
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)


temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
tab <- TabixFile(zipped, idx, yieldSize = 8000)
open(tab)
while(nrow(vcf <- readVcf(tab, genome = opt$genome))) {
    exptData(vcf)$header <- newVcfHeader

    oldwd <- getwd()

    if (sum(rowData(vcf)$FILTER == "PASS") > 0) {
        # convert vcf to chasm input
        al <- sapply(rowData(vcf)$ALT, length)
        X <- data.frame(id = 1:nrow(vcf), seq = as.character(seqnames(rowData(vcf))), zero = as.integer(start(rowData(vcf)) - 1), one = as.integer(start(rowData(vcf))), strand = rep('+', nrow(vcf)), ref = as.character(rowData(vcf)$REF), filt = rowData(vcf)$FILTER, stringsAsFactors = F)
        re <- rep.int(X$id, times = al)
        X <- X[re, ,drop = F]
        X$alt <- as.character(unlist(rowData(vcf)$ALT))
        X <- subset(X, sapply(alt, nchar) == 1 & sapply(ref, nchar) == 1 & filt == "PASS")
        X <- X[,-which(colnames(X) == 'filt')]
        X$seq <- sub('^', 'chr', X$seq)

        setwd(opt$chasmDir)
        tmp <- tempfile()
        write.table(X, file = tmp, quote = F, sep = '\t', row.names = F, col.names = F)
        for (cl in classifiers) {
            cmd <- paste("CHASMDIR=", opt$chasmDir, " ", opt$python, ' ./RunChasm ', cl, ' ', tmp, ' -g' ,sep = '')
            #cat(cmd)
            system(cmd, ignore.stdout = T)
            results <- read.table(file = paste(tmp, cl, '.output', sep = ''), sep = '\t', header = T, as.is = T)
            if (nrow(results) > 1) {
                infoprime <- info(vcf)
                infoprime[as.integer(as.character(results$MutationID)), paste(cl, "chasm_mut", sep = "_")] <- as.character(results$Mutation)
                infoprime[as.integer(as.character(results$MutationID)), paste(cl, "chasm_score", sep = "_")] <- results$CHASM
                infoprime[as.integer(as.character(results$MutationID)), paste(cl, "chasm_pval", sep = "_")] <- results$PValue
                if (nrow(results) > 10) {
                    infoprime[as.integer(as.character(results$MutationID)), paste(cl, "chasm_fdr", sep = "_")] <- results$BHFDR
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
