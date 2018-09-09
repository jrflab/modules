#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

optList <- list(make_option("--outFile", default = NULL, help = "output file"));
parser <- OptionParser(usage = "%prog [options] mutation_file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input mutation table files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;
if (length(files)>1) { cat ("Only converting first file. Everything else is ignored!!!!!\n")}

dat <- read.delim(files[1], as.is=T)

res <- lapply(unique(dat$cluster_id), function(cl) {
	smalldat <- subset(dat, cluster_id==cl)
	sample_levels=sort(unique(smalldat$sample_id))
	smalldat$sample_id <- factor(smalldat$sample_id, levels=sample_levels)
	data.frame(
		sample_id = sample_levels,
		cluster_id=cl,
		size=as.numeric(table(smalldat$sample_id )),
		mean=tapply(smalldat$cellular_prevalence, smalldat$sample_id , mean),
		std=tapply(smalldat$cellular_prevalence, smalldat$sample_id , sd))
	})
write.table(do.call("rbind", res), file=opt$outFile, sep="\t", row.names=F, na="", quote=F)
