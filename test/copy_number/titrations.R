#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--tumor_normal", default = NA, type = 'character', help = "tumor and normal samples"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_sample = strsplit(x=opt$tumor_normal, split="_", fixed=TRUE)[[1]][1]
normal_sample = strsplit(x=opt$tumor_normal, split="_", fixed=TRUE)[[1]][2]
load(paste0("facets/cncf/", opt$tumor_normal, ".Rdata"))
purity = c(.1, .2, .3, .4, .5, .75, .9)
purity = purity[purity < fit$purity]
for (i in 1:length(purity)) {
	system(paste0("samtools view -s 0", purity[i], " -b bam/", tumor_sample, ".bam > titration/", tumor_sample, purity[i], ".bam"))
	system(paste0("samtools merge titration/", tumor_sample, purity[i], ".bam bam/", normal_sample[i], ".bam titration/", tumor_sample, "-", purity[i], ".bam"))
	file.remove(paste0("titration/", tumor_sample, purity[i], ".bam"))
	cat(sessionInfo()$R.version$version.string, file=paste0("titrations/", opt$tumor_normal, ".timestamp"))
}
