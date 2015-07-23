#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

suppressPackageStartupMessages(library("optparse"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--sampleName", action = "store_true", default = F, help = "add samplename column [default %default]"),
                make_option("--normalLast", action = "store_true", default = F, help = "normal sample last [default %default]"),
                make_option("--tumorNormal", action = "store_true", default = F, help = "add tumor-normal samplename column [default %default]"))

parser <- OptionParser(usage = "%prog [options] vcf.file", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;
files <- arguments$args;


Data <- list();
for (f in files) {
    X <- read.table(file = f, sep = '\t', as.is = T, comment.char = '', quote = '');
    if (nrow(X) <= 1) {
        next
    }
    h <- X[1,];
    h <- sub('#', '', h)
    colnames(X) <- h
    X <- X[-1, ]

    if (opt$sampleName) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        X[,"SAMPLE"] <- sname
        h <- c(h, "SAMPLE")
        h <- sub(paste(sname, '\\.', sep = ''), 'SAMPLE.', h)
        colnames(X) <- h
        Data[[sname]] <- X
    }
    if (opt$tumorNormal) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        tumor <- sub('_.*', '', sname)
        normal <- sub('.*_', '', sname)
        X[,"TUMOR_SAMPLE"] <- tumor
        h <- c(h, "TUMOR_SAMPLE")
        X[,"NORMAL_SAMPLE"] <- normal
        h <- c(h, "NORMAL_SAMPLE")

        h <- sub(paste(tumor, '\\.', sep = ''), 'TUMOR.', h)
        h <- sub(paste(normal, '\\.', sep = ''), 'NORMAL.', h)
        colnames(X) <- h
        Data[[sname]] <- X
    }
    if (opt$normalLast) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        normal <- sub('.*_', '', sname)
        normFields <- grep(paste(normal, ".", sep = ""), h, fixed = T, value = T)
        fields <- sub('.*\\.', '', normFields)
        x <- grep(paste("\\.", fields[1], "$", sep = ""), h, perl = T, value = T)
        samples <- sub('\\..*', '', x)
        tumorSamples <- samples[-which(samples == normal)]
        for (i in 1:length(tumorSamples)) {
            tumor <- tumorSamples[i]
            hh <- h
            XX <- X
            for (otherTumor in tumorSamples[-i]) {
                x <- grep(paste("^", otherTumor, "\\.", sep = ""), hh, perl = T)
                hh <- hh[-x]
                XX <- XX[, -x]
            }
            XX[,"TUMOR_SAMPLE"] <- tumor
            hh <- c(hh, "TUMOR_SAMPLE")
            XX[,"NORMAL_SAMPLE"] <- normal
            hh <- c(hh, "NORMAL_SAMPLE")
            hh <- sub(paste(tumor, '\\.', sep = ''), 'TUMOR.', hh)
            hh <- sub(paste(normal, '\\.', sep = ''), 'NORMAL.', hh)
            colnames(XX) <- hh
            Data[[tumor]] <- XX
        }
    }   
}
if (length(Data) == 0) {
    # print empty table in case there are no mutations
    f <- files[1]
    X <- read.table(file = f, sep = '\t', as.is = T, comment.char = '', quote = '')
    h <- X[1,]
    h <- sub('#', '', h)
    colnames(X) <- h
    X <- X[-1, ]
    if (opt$sampleName) {
        X[,"SAMPLE"] <- character(0)
        h <- c("SAMPLE", h)
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        h <- sub(paste(sname, '\\.', sep = ''), 'SAMPLE.', h)
        colnames(X) <- h
    }
    if (opt$tumorNormal) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        tumor <- sub('_.*', '', sname)
        normal <- sub('.*_', '', sname)
        X[,"TUMOR_SAMPLE"] <- character(0)
        X[,"NORMAL_SAMPLE"] <- character(0)
        h <- c("TUMOR_SAMPLE", "NORMAL_SAMPLE", h)

        h <- sub(paste(tumor, '\\.', sep = ''), 'TUMOR.', h)
        h <- sub(paste(normal, '\\.', sep = ''), 'NORMAL.', h)
        colnames(X) <- h
    }
    write.table(X, sep = '\t', row.names = F, quote = F)

    quit(save = 'no', status = 0)
}

fields <- unique(unlist(lapply(Data, colnames)))

for (f in names(Data)) {
    miss <- setdiff(fields, colnames(Data[[f]]));
    Data[[f]][,miss] <- NA;
    Data[[f]] <- Data[[f]][, fields];
}
table.merged <- do.call(rbind, Data);
rownames(table.merged) <- NULL

if (opt$sampleName) {
    x <- which(colnames(table.merged) == "SAMPLE")
    y <- which(colnames(table.merged) != "SAMPLE")
    table.merged <- table.merged[, c(x,y)]
}
if (opt$tumorNormal || opt$normalLast) {
    xx <- colnames(table.merged) == "TUMOR_SAMPLE" | colnames(table.merged) == "NORMAL_SAMPLE"
    x <- which(xx)
    y <- which(!xx)
    table.merged <- table.merged[, c(x,y)]
}

write.table(table.merged, sep = '\t', row.names = F, quote = F)
