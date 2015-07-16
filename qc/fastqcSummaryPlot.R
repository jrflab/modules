#!/usr/bin/env Rscript

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("ggplot2"));

optList <- list(
                make_option("--outPrefix", default = 'summary', type = "character", action = "store", help ="Output file prefix (default %default)"),
                make_option("--width", default = 1000, action = "store", help ="width of heatmap image (default %default)"),
                make_option("--height", default = 1000, action = "store", help ="height of heatmap image (default %default)"))
parser <- OptionParser(usage = "%prog [options] summary file(s)", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

summ <- arguments$args %>%
    llply(read.table, sep = '\t', stringsAsFactors = F) %>%
    bind_rows %>%
    setNames(c("Status", "Metric", "Sample")) %>%
    mutate(Sample = str_replace(Sample, '\\.bam', '')) %>%
    mutate_each(funs(as.factor))

fn <- str_c(opt$outPrefix, '.png')
png(fn, width = opt$width + 10 * nlevels(summ$Sample), height = opt$height, type = 'cairo-png');
p <- summ %>% ggplot(aes(x = Sample, y = Metric)) +
    geom_tile(aes(fill = Status)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off();

fn <- str_c(opt$outPrefix, '.txt')
summ %>% spread(Metric, Status) %>% write.table(file = fn, sep = '\t', quote = F, row.names = F)
