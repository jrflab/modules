#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--sample_set", default = NULL, help = "sample set name"),
			   make_option("--normal_samples", default = NULL, help = "normal sample names"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_samples = unlist(strsplit(opt$sample_set, split="_", fixed=TRUE))
normal_sample = unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE))
normal_sample = tumor_samples[tumor_samples %in% normal_sample]
tumor_samples = tumor_samples[!(tumor_samples %in% normal_sample)]

cat("num_iters: 10000\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = FALSE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("base_measure_params:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  alpha: 1\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  beta: 1\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("concentration:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  value: 1.0\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  prior:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("    shape: 1.0\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("    rate: 0.001\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("density: pyclone_beta_binomial\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("beta_binomial_precision_params:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  value: 1000\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  prior:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("    shape: 1.0\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("    rate: 0.0001\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("  proposal:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("    precision: 0.5\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat(paste0("working_dir: pyclone/",opt$sample_set, "\n"), file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("trace_dir: trace", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("init_method: connected\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
cat("samples:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)

for (i in 1:length(tumor_samples)) {
	if (i!=1) {
		cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	}
	cat(paste0("  ", tumor_samples[i], ":\n"), file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	cat(paste0("    mutations_file: ", tumor_samples[i], ".yaml\n"), file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	cat("    tumour_content:\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	load(paste0("ascat/ascat/", tumor_samples[i], "_", normal_sample, ".RData"))
	cat(paste0("      value: ", ifelse(is.na(purity), 1.0, signif(purity, 2)),"\n"), file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	cat("    error_rate: 0.01", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	if (i!=length(tumor_samples)) {
		cat("\n", file=paste0("pyclone/", opt$sample_set, "/config.yaml"), append = TRUE)
	}
}

for (i in 1:length(tumor_samples)) {
	system(paste0("source ~/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate ~/share/usr/anaconda-envs/PyClone-0.13.1 && PyClone build_mutations_file --in_file pyclone/",  opt$sample_set, "/", tumor_samples[i], ".tsv --out_file pyclone/", opt$sample_set, "/", tumor_samples[i], ".yaml  --prior parental_copy_number"))
}
