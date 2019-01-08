#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("VariantAnnotation"))

#options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
#options(error = recover)
options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'b37', help = "genome build [default %default]"),
                make_option("--tumor", default = NULL, help = "tumor sample"),
                make_option("--outFile", default = NULL, help = "output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.files", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) <= 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

vcfs <- list()
for (f in files) {
    vcf <- readVcf(f, genome = opt$genome)
    tum <- ifelse(opt$tumor %in% colnames(geno(vcf)$GT), opt$tumor, "TUMOR")
    gt <- geno(vcf)$GT[, tum]
    vcf <- vcf[gt != "./." & gt != "0/0" & gt != "0", ]
    vcfs <- append(vcfs, vcf)
}

all <- do.call('rbind', lapply(vcfs, function(x) as.data.frame(subset(rowRanges(x), FILTER == "PASS"))))
all <- all[, c("seqnames", "start", "end")]
cnt <- ddply(all, .(seqnames, start, end), nrow)

cnt <- subset(cnt, V1 > 1)
write(paste(cnt[,1], ":", cnt[,2], "-", cnt[,3], sep = ''), file = opt$outFile)
#write.table(cnt[,c(1:3)], file = opt$outFile, sep = '\t', quote = F, row.names = F, col.names = F)
