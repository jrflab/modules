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
t_fraction = c(".10", ".20", ".30", ".40", ".50", ".60", ".70", ".80", ".90")
n_fraction = rev(t_fraction)
for (i in 1:length(t_fraction)) {
	system(paste0("samtools view -s 0", t_fraction[i], " -b bam/", tumor_sample, ".bam > titrations/bam/", tumor_sample, t_fraction[i], ".bam"))
	system(paste0("samtools view -s 0", n_fraction[i], " -b bam/", normal_sample, ".bam > titrations/bam/", normal_sample, n_fraction[i], ".bam"))
	system(paste0("samtools merge  titrations/bam/", tumor_sample, t_fraction[i], n_fraction[i], ".bam titrations/bam/", tumor_sample, t_fraction[i], ".bam titrations/bam/", normal_sample, n_fraction[i],".bam"))
	system(paste0("samtools index titrations/bam/", tumor_sample, t_fraction[i], n_fraction[i], ".bam"))
	file.copy(from=paste0("titrations/bam/", tumor_sample, t_fraction[i], n_fraction[i], ".bam.bai"), to=paste0("titrations/bam/", tumor_sample, t_fraction[i], n_fraction[i], ".bai"))
	file.remove(paste0("titrations/bam/", tumor_sample, t_fraction[i], ".bam"))
	file.remove(paste0("titrations/bam/", normal_sample, n_fraction[i], ".bam"))
}
cat(sessionInfo()$R.version$version.string, file=paste0("titrations/bam/", opt$tumor_normal, ".timestamp"))
