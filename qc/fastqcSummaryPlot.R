#!/usr/bin/env Rscript

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

suppressPackageStartupMessages(library("optparse"));

optList <- list(
                make_option("--outFile", default = 'summary.png', type = "character", action = "store", help ="Output file (default %default)"),
                make_option("--width", default = 1000, action = "store", help ="width of heatmap image (default %default)"),
                make_option("--height", default = 1000, action = "store", help ="height of heatmap image (default %default)"))
parser <- OptionParser(usage = "%prog [options] summaryFIle", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

summFile <- arguments$args[1];
summ <- as.matrix(read.table(summFile, row.names = 1, sep = '\t', header = T, as.is = T));
summ[summ == "PASS"] <- 1;
summ[summ == "FAIL"] <- 2;
summ[summ == "WARN"] <- 3;
class(summ) <- 'numeric';

png(opt$outFile, width = opt$width, height = opt$height, type = 'cairo-png');
heatmap(t(summ), col = c("black", "red", "yellow"), margins = c(10, 20), scale = 'none');
legend('topleft', legend = c("pass", "fail", "warn"), fill = c('black', 'red', 'yellow'));
null <- dev.off();



