#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))

optList = list(make_option("--sample_name", default = NULL, help = "tumor normal sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_sample = unlist(strsplit(opt$sample_name, split="_", fixed=TRUE))[1]
normal_sample = unlist(strsplit(opt$sample_name, split="_", fixed=TRUE))[2]

mutation_summary = read_tsv(file="summary/tsv/mutation_summary.tsv", col_types = cols(.default = col_character())) %>%
				   type_convert() %>%
				   filter(TUMOR_SAMPLE==tumor_sample) %>%
 				   filter(NORMAL_SAMPLE==normal_sample) %>%
 				   filter(grepl("mutect", variantCaller, fixed=TRUE)) %>%
 				   filter(NORMAL_MAF==0) %>%
 				   filter(TUMOR_MAF>=.05) %>%
 				   filter(TUMOR_DP<=500) %>%
 				   filter(TUMOR_DP>=20) %>%
 				   filter(NORMAL_DP<=500) %>%
 				   filter(NORMAL_DP>=10) %>%
 				   mutate(CHROM = as.numeric(ifelse(CHROM=="X", 23, CHROM))) %>%
 				   mutate(CHROM = as.numeric(ifelse(CHROM=="Y", 24, CHROM))) %>%
 				   filter(CHROM<=22) %>%
 				   mutate(UUID = paste0(CHROM, ":", POS, "_", REF, "_", ALT))

load(paste0("facets/cncf/", opt$sample_name, ".Rdata"))
qt = q1 = rep(NA, nrow(mutation_summary))
for (i in 1:nrow(mutation_summary)) {
	x = mutation_summary$CHROM[i]
	y = mutation_summary$POS[i]
	indx = which(fit$cncf[,"chrom"]==x & (fit$cncf[,"start"]<=y & fit$cncf[,"end"]>=y))
	if (length(indx)!=0) {
		qt[i] = fit$cncf[indx,"tcn.em"]
		q1[i] = fit$cncf[indx,"lcn.em"]
	}
}
fsq = as.numeric(mutation_summary$TUMOR_MAF)
n = as.numeric(mutation_summary$TUMOR_DP)
mutation_id = as.character(mutation_summary$UUID)
var_counts = round(fsq*n)
ref_counts = round((1-fsq)*n)
normal_cn = rep(2, nrow(mutation_summary))
minor_cn = q1
major_cn = qt-q1
sample_summary = data.frame(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)
index = apply(sample_summary, 1, function(x) {any(is.na(x))})
sample_summary = sample_summary[!index,,drop=FALSE]
index = sample_summary[,"major_cn"]==0
sample_summary = sample_summary[!index,,drop=FALSE]
write.table(sample_summary, paste0("pyclone/", opt$sample_name, "/", tumor_sample,".tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)

cat("num_iters: 100000\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = FALSE)
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

cat(paste0("  ", tumor_sample, ":\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat(paste0("    mutations_file: ", tumor_sample, ".yaml\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    tumour_content:\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat(paste0("      value: ", ifelse(is.na(fit$purity), 1.0, signif(fit$purity, 2)),"\n"), file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("\n", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
cat("    error_rate: 0.01", file=paste0("pyclone/", opt$sample_name, "/config.yaml"), append = TRUE)
system(paste0("source ~/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate ~/share/usr/anaconda-envs/PyClone-0.13.1 && PyClone build_mutations_file --in_file pyclone/",  opt$sample_name, "/", tumor_sample, ".tsv --out_file pyclone/", opt$sample_name, "/", tumor_sample, ".yaml  --prior parental_copy_number"))
