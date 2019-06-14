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

in_file = list(
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/", tumor_sample,".tsv"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/pyclone.tsv")
)
out_file = list(
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_std_by_cid.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_ccf_by_cid.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_std_by_cn.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_ccf_by_cn.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/summary.tsv")
)

mutation_summary = read_tsv(file=in_file[[1]], col_types = cols(.default = col_character())) %>%
				   type_convert()

pyclone_summary = read_tsv(file=in_file[[2]], col_types = cols(.default = col_character()), col_names = c("mutation_id", "ccf", "std", "cluster_id")) %>%
				  type_convert() %>%
				  mutate(cluster_id = factor(cluster_id)) %>%
				  mutate(ccf = as.numeric(ccf)) %>%
				  mutate(std = as.numeric(std)) %>%
				  slice(-1)
				  
mutation_summary = full_join(mutation_summary, pyclone_summary, by="mutation_id")

plot.0 =  ggplot(mutation_summary, aes(x=std, fill=cluster_id)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x=expression(sigma), y="Frequency\n") +
		  guides(fill=guide_legend(title=c("Cluster")))
pdf(file=out_file[[1]], width=6, height=6)
print(plot.0)
dev.off()
		 
plot.0 =  ggplot(mutation_summary, aes(x=ccf, fill=cluster_id)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="CCF", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  guides(fill=guide_legend(title=c("Cluster")))
pdf(file=out_file[[2]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x=std, fill=major_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x=expression(sigma), y="Frequency\n") +
		  coord_cartesian(xlim=c(0,100)) +
		  guides(fill=guide_legend(title=c("Copy number")))

pdf(file=out_file[[3]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x=ccf, fill=major_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="CCF", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  guides(fill=guide_legend(title=c("Copy number")))

pdf(file=out_file[[4]], width=6, height=6)
print(plot.0)
dev.off()

write_tsv(x=mutation_summary, path=out_file[[5]])
