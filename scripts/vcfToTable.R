#!/usr/bin/env Rscript
# converts a vcf file to a tab-delimited table that can be more easily parsed

suppressPackageStartupMessages(library("optparse"));

optList <- list(
                make_option("--outFile", default = NULL, help = "Output file [default STDOUT]"),
                make_option("--cbindSamples", action = "store_true", default = F, help = "Col-concatenate sample chromosome position calls instead of row concatenating calls [default %default]"),
                #make_option("--showFiltered", action = "store_true", default = F, help = "Show filtered records [default %default]"),
                #make_option("--includeHomRef", action = "store_true", default = F, help = "Include homozygous reference entries [default %default]"),
                #make_option("--includeFiltered", action = "store_true", default = F, help = "Include filtered entries [default %default]"),
                #make_option("--onlyNovel", action = "store_true", default = F, help = "Only include novel entries [default %default]"),
                make_option("--fields", default = NULL, help = "File containing fields to output [default ALL]"));

parser <- OptionParser(usage = "%prog [options] vcf.file", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need input VCF file\n");
    print_help(parser);
    stop();
} else {
    vcf.file <- arguments$args[1];
}

if (!is.null(opt$fields)) {
    userFields <- scan(opt$fields, what = 'character', quiet = T);
}

f <- file(vcf.file, "r");
header <- readLines(f, n = 1000);
header <- grep('^#', header, perl = T, value = T);
header <- header[length(header)];
header <- sub('^#', '', header, perl = T);
header <- strsplit(header, "\t")[[1]];
close(f);
vcf <- read.table(vcf.file, sep = '\t', comment.char = "#", quote = '', as.is = T);
colnames(vcf) <- make.names(header);

if (ncol(vcf) >= 9 && colnames(vcf)[9] == "FORMAT") {
    samples <- colnames(vcf)[10:ncol(vcf)]

    # gather fields
    ff <- unique(vcf$FORMAT);
    if (length(ff) > 1) {
        cat("> 1 unique format field not supported");
        stop();
    }
    formatFields <- unlist(strsplit(ff, ':'));
} else {
    samples <- NULL
    formatFields <- c()
}

infoList <- sapply(sapply(vcf$INFO, strsplit, ";"), function(x) {
               sp <- strsplit(x, "=");
               vals <- sapply(sp, function(x) x[2]);
               names(vals) <- sapply(sp, function(x) x[1]);
               vals;
                });
names(infoList) <- NULL;

infoFields <- unique(unlist(c(sapply(infoList, names))));
info <- do.call(rbind, lapply(lapply(infoList, unlist), "[", infoFields));
colnames(info) <- infoFields;

vcfFields <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER");
fields <- c("SAMPLE", vcfFields, infoFields, formatFields);

if (!is.null(opt$fields)) {
    if (any(!userFields %in% fields)) {
        cat("Missing field in vcf file");
        stop();
    }
    fields <- unique(c("SAMPLE", userFields))
}

if (!is.null(samples)) {
    Data <- list()
    if (opt$cbindSamples) {
        # one row per position, multiple fields per sample
        sampleGT <- list()
        for (sp in samples) {
            sampleGT[[sp]] <- do.call(rbind, sapply(vcf[,sp], strsplit, ':'));
            rownames(sampleGT[[sp]]) <- NULL;
            colnames(sampleGT[[sp]]) <- formatFields;
            colnames(sampleGT[[sp]]) <- paste(sp, colnames(sampleGT[[sp]]), sep = '.')
        }
        Data <- cbind(vcf[, vcfFields], do.call(cbind, sampleGT), info, stringsAsFactors = F);
    } else {
        # row concat sample calls
        for (sp in samples) {
            sampleGT <- sapply(vcf[, sp], strsplit, ':');
            sampleGT <- do.call(rbind, sampleGT);
            rownames(sampleGT) <- NULL;
            colnames(sampleGT) <- formatFields;
            colnames(sampleGT) <- paste("Sample", colnames(sampleGT), sep = '')
            sampleInfo <- cbind(SAMPLE = rep(sp, nrow(sampleGT)), vcf[, vcfFields], sampleGT, info, stringsAsFactors = F);
            #sampleInfo <- sampleInfo[sampleInfo[, "SampleGT"] != "./.", , drop = F];
            #if (!opt$includeFiltered) {
                #sampleInfo <- sampleInfo[sampleInfo[, "FILTER"] == "PASS", , drop = F];
            #}
            #if (!opt$includeHomRef) { 
                #sampleInfo <- sampleInfo[sampleInfo[, "SampleGT"] != "0/0", , drop = F];
            #}
            Data[[sp]] <- sampleInfo;
        }
        Data <- do.call(rbind, Data);
    }
    Data <- Data[, !apply(Data, 2, function(x) all(is.na(x)))];
    Data <- Data[, !duplicated(colnames(Data))];
} else {
    Data <- cbind(vcf[, vcfFields], info, stringsAsFactors = F)
}

if (is.null(opt$outFile)) {
    out <- "";
} else {
    out <- opt$outFile;
}

write.table(Data, file = out, quote = F, sep = "\t", row.names = F)

