#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("VariantAnnotation"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--maData", default = NULL, help = "mutation assessor  R data"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]"))
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$maData)) {
    cat("Need Mutation assessor Rdata file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}
vcfFile <- arguments$args[1];

load(opt$maData)
#maFiles <- paste("MA.chr", c(sapply(1:22, function (x) sprintf("%02d", x)), "X", "Y"), ".txt", sep = "")
#names(maFiles) <- c(1:22, 'X', 'Y')
## load mutation assessor database
#mutacc <- list()
#oldwd <- getwd()
#setwd(opt$maDir)
#for (chr in names(maFiles)) {
#    fn <- maFiles[chr]
#    cat("Reading", fn, "\n")
#    mutacc[[chr]] <- read.table(fn, sep = '\t', header = T, quote = '', comment.char = '', as.is = T)
#}
#setwd(oldwd)

vcf <- readVcf(vcfFile, genome = 'hg19')

hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
rownames(hinfoprime) <- rownames(info(header(vcf)))
hinfoprime <- rbind(hinfoprime, mut_ass_refgen_variant = c("A", "String", "RefGenome variant"))
hinfoprime <- rbind(hinfoprime, mut_ass_gene = c("A", "String", "Mutation assessor Gene"))
hinfoprime <- rbind(hinfoprime, mut_ass_uniprot = c("A", "String", "Mutation assessor Uniprot"))
hinfoprime <- rbind(hinfoprime, mut_ass_info = c("A", "String", "Mutation assessor Info"))
hinfoprime <- rbind(hinfoprime, mut_ass_uniprot_var = c("A", "String", "Mutation assessor Uniprot variant"))
hinfoprime <- rbind(hinfoprime, mut_ass_impact = c("A", "String", "Mutation assessor functional impact"))
hinfoprime <- rbind(hinfoprime, mut_ass_score = c("A", "Float", "Mutation assessor score"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(metadata(vcf)$header)
hlist$INFO <- hinfoprime
metadata(vcf)$header <- new("VCFHeader", samples = header(vcf)@samples, header = hlist)


n <- sapply(rowRanges(vcf)$ALT, length)
chr <- factor(rep(as.character(seqnames(rowRanges(vcf))), n))
x <- data.frame(rowId = rep(1:nrow(vcf), n), mutAssId = paste(rep('hg19', sum(n)), rep(as.character(seqnames(rowRanges(vcf))), n), rep(start(rowRanges(vcf)), n), rep(rowRanges(vcf)$REF, n), unlist(rowRanges(vcf)$ALT), sep = ','), stringsAsFactors = F)
x <- split(x, chr)

for (ch in names(x)) {
    x[[ch]] <- merge(x[[ch]], mutass[[ch]], by.x = 'mutass', by.y = 'Mutation')
}
mutAssCols <- c('mut_ass_refgen_variant', 'mut_ass_gene', 'mut_ass_uniprot', 'mut_ass_info', 'mut_ass_uniprot_var', 'mut_ass_impact', 'mut_ass_score')
results <- as.data.frame(rbindlist(x))
colnames(results) <- c('mutAssId', 'rowId', mutAssCols)

infoprime <- info(vcf)
infoprime[, mutAssCols] <- as.character(rep(NA, nrow(vcf)))
infoprime[, "mut_ass_score"] <- as.numeric(rep(NA, nrow(vcf)))
info(vcf) <- infoprime

if (!is.null(results) && nrow(results) > 0) {
    cat("Merging mutAss results ... ")
    for (co in mutAssCols) {
        info(vcf)[results$rowId, co] <- results[, co]
    }
} else {
    cat("No results from fathmm\n")
}

#fix sample genotype order

x <- which(names(geno(vcf)) == "GT")
ord <- c(x, (1:length(geno(vcf)))[-x])
geno(vcf) <- geno(vcf)[ord]

cat("Writing vcf to", opt$outFile, "... ")
setwd(oldwd)
writeVcf(vcf, opt$outFile)
cat("done\n")

