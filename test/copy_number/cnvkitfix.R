#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--tumor", default = NA, type = 'character', help = "tumor sample"),
				 make_option("--normals", default = NA, type = 'character', help = "normal samples"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_sample = opt$tumor
normal_samples = unlist(strsplit(x=opt$normals, " ", fixed=TRUE))

file.copy(from=paste0("cnvkit/cnn/tumor/", tumor_sample, ".antitargetcoverage.cnn"), to=paste0("cnvkit/cnn/tumor/", tumor_sample, ".A.antitargetcoverage.cnn"))
for (i in 1:length(normal_samples)) {
	system(paste0("cnvkit.py fix cnvkit/cnn/tumor/", tumor_sample, ".A.targetcoverage.cnn  cnvkit/cnn/tumor/", tumor_sample, ".A.antitargetcoverage.cnn cnvkit/reference/", normal_samples[i], ".A.cnr -o cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".A.cnr"))
}
file.remove(paste0("cnvkit/cnn/tumor/", tumor_sample, ".A.antitargetcoverage.cnn"))

file.copy(from=paste0("cnvkit/cnn/tumor/", tumor_sample, ".antitargetcoverage.cnn"), to=paste0("cnvkit/cnn/tumor/", tumor_sample, ".B.antitargetcoverage.cnn"))
for (i in 1:length(normal_samples)) {
	system(paste0("cnvkit.py fix cnvkit/cnn/tumor/", tumor_sample, ".B.targetcoverage.cnn  cnvkit/cnn/tumor/", tumor_sample, ".B.antitargetcoverage.cnn cnvkit/reference/", normal_samples[i], ".B.cnr -o cnvkit/cnr/", tumor_sample, "_", normal_samples[i], ".B.cnr"))
}
file.remove(paste0("cnvkit/cnn/tumor/", tumor_sample, ".B.antitargetcoverage.cnn"))
