#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))

optList = list(
			   make_option("--sample_set", default = NULL, help = "sample set name"),
			   make_option("--normal_samples", default = NULL, help = "normal sample names"),
			   make_option("--min_depth", default = NA, help = "minimum depth to consider")
			   )
			   
parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

tumor_samples = unlist(strsplit(opt$sample_set, split="_", fixed=TRUE))
normal_sample = unlist(strsplit(opt$normal_samples, split=" ", fixed=TRUE))
normal_sample = tumor_samples[tumor_samples %in% normal_sample]
tumor_samples = tumor_samples[!(tumor_samples %in% normal_sample)]
min_depth = ifelse(is.na(opt$min_depth) | is.null(opt$min_depth) | opt$min_depth=="" | opt$min_depth==" ", 50, opt$min_depth)

mutation_summary = read_tsv(file=paste0("sufam/", opt$sample_set, ".tsv"))
index = apply(mutation_summary[,paste0("DP_", tumor_samples)], 1, function(x) {sum(x>=min_depth)})==length(tumor_samples)
mutation_summary = mutation_summary[index,,drop=FALSE]
pyclone_summary = read_tsv(file=paste0("pyclone/", opt$sample_set, "/report/pyclone.tsv"), col_types = cols(.default = col_character())) %>%
				  type_convert() %>%
				  bind_cols(mutation_summary)
				  
clusters = table(pyclone_summary$cluster_id)
if (any(clusters==1)) {
	pyclone_summary = pyclone_summary %>%
					  filter(!(cluster_id %in% names(clusters)[clusters==1]))
}

pdf(file=paste0("pyclone/", opt$sample_set, "/report/pyclone.pdf"), width=7, height=6)
for (i in 1:(length(tumor_samples)-1)) {
	for (j in (i+1):length(tumor_samples)) {
		x = pyclone_summary[,tumor_samples[i]] %>%
			.[[1]]
		y = pyclone_summary[,tumor_samples[j]] %>%
			.[[1]]
		z = pyclone_summary %>%
			.[["cluster_id"]]
		c_x = pyclone_summary %>%
		  	  .[[paste0("CALL_", tumor_samples[i])]]
		c_y = pyclone_summary %>%
		  	  .[[paste0("CALL_", tumor_samples[j])]]
		x[c_x==0] = 0
		y[c_y==0] = 0
		tmp.0 = data_frame(x=x, y=y, z=factor(z, ordered=TRUE))
		plot.0 =  ggplot(tmp.0, aes(x=x, y=y, fill=z, color=z)) +
				  geom_point(alpha = .8, size=2) +
				  geom_contour() +
				  theme_classic() +
				  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
				  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
				  labs(x=paste0("\n",tumor_samples[i],"\n"), y=paste0("\n",tumor_samples[j],"\n")) +
				  guides(color=guide_legend(title=c("Cluster"))) +
				  guides(fill=FALSE)
		print(plot.0)
	}
}
dev.off()
