#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

optList = list(make_option("--sample_name", default = NULL, help = "sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

file_names = dir(paste0("ascat/ascat"), pattern=".RData", full.names=FALSE)
file_names = gsub(paste0("_", opt$sample_name), "", x=gsub(".RData", "", x=file_names[grep(opt$sample_name, file_names)], fixed=TRUE), fixed=TRUE)

cat("num_iters: 10000\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = FALSE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("base_measure_params:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  alpha: 1\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  beta: 1\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("concentration:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  value: 1.0\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  prior:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    shape: 1.0\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    rate: 0.001\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("density: pyclone_beta_binomial\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("beta_binomial_precision_params:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  value: 1000\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  prior:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    shape: 1.0\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    rate: 0.0001\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("  proposal:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    precision: 0.5\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat(paste0("working_dir: pyclone/",opt$sample_name, "\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("trace_dir: trace", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("init_method: connected\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("samples:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)

for (i in 1:length(file_names)) {
	if (i!=1) {
		cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	}
	cat(paste0("  ", file_names[i], ":\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	cat(paste0("    mutations_file: ", file_names[i], ".yaml\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	cat("    tumour_content:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	load(paste0("ascat/ascat/", file_names[i], "_", opt$sample_name, ".RData"))
	cat(paste0("      value: ", ifelse(is.na(purity), 1.0, signif(purity, 2)),"\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	cat("    error_rate: 0.01", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	if (i!=length(file_names)) {
		cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
	}
}

