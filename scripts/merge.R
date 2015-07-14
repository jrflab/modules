#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

optList <- list(
                make_option(c("-a", "--all"), action = "store_true", default = F, help = "shorthand for all.x and all.y"),
                make_option(c("-b", "--byCommon"), action = "store_true", default = F, help ="join on common columns"),
                make_option(c("-y", "--byY"), default = NULL, help ="common column in y"),
                make_option(c("-x", "--byX"), default = NULL, help ="common column in x"),
                make_option("--byColY", default = NULL, help ="common column y"),
                make_option("--byColX", default = NULL, help ="common column x"),
                make_option(c("-X", "--allX"), action = "store_true", default = F, help = "outer join on x"),
                make_option(c("-Y", "--allY"), action = "store_true", default = F, help = "outer join on y"),
                make_option(c("-i", "--incomparables"), default = NA, help = "Values which cannot be matched [default %default]"),
                make_option(c("-s", "--sep"), default = "\t", help = "field separator [default \\t]"),
                make_option(c("-H", "--header"), action = "store_true", default = F, help = "header"));


parser <- OptionParser(usage = "%prog [options] file1 file2", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need two files to merge\n\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

for (file in files) {
    if (file.access(file) == -1) {
        stop(sprintf("Specified file ( %s ) does not exist", file));
    }
}

X <- read.table(files[1], sep = opt$sep, as.is = T, header = opt$header, quote = '');
Y <- read.table(files[2], sep = opt$sep, as.is = T, header = opt$header, quote = '');

if (opt$byCommon) {
    opt$by <- intersect(names(X), names(Y));
}

if (is.null(opt$byY)) {
    opt$byY <- opt$by;
}

if (is.null(opt$byX)) {
    opt$byX <- opt$by;
}

if (!is.null(opt$byColY)) {
    opt$byY <- as.integer(opt$byColY);
}

if (!is.null(opt$byColX)) {
    opt$byX <- as.integer(opt$byColX);
}
m <- merge(X, Y, by = opt$by, by.x = opt$byX, by.y = opt$byY, all = opt$all, all.x = opt$allX, all.y = opt$allY, incomparables = opt$incomparables);

write.table(m, sep = opt$sep, quote = F, row.names = F);
