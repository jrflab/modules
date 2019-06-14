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
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_vaf_by_cn.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_depth_by_cn.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/scatter_vaf_depth_by_cn.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/summary.tsv"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/clusters.tsv")
)

mutation_summary = read_tsv(file=in_file[[1]], col_types = cols(.default = col_character())) %>%
				   type_convert() %>%
				   mutate(total_cn = factor(minor_cn+major_cn)) %>%
				   mutate(DP = var_counts+ref_counts) %>%
				   mutate(VAF = 100*var_counts/(var_counts+ref_counts))

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
		  labs(x="\nCCF\n", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  guides(fill=guide_legend(title=c("Cluster")))
pdf(file=out_file[[2]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x=std, fill=total_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x=expression(sigma), y="Frequency\n") +
		  guides(fill=guide_legend(title=c("Copy number")))

pdf(file=out_file[[3]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x=ccf, fill=total_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nCCF\n", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  guides(fill=guide_legend(title=c("Copy number")))

pdf(file=out_file[[4]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x = VAF, fill=total_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nVAF(%)\n", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,100)) +
		  guides(fill=guide_legend(title=c("Copy number")))
		  
pdf(file=out_file[[5]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x = DP, fill=total_cn)) +
		  geom_histogram(alpha = .8) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nDP\n", y="Frequency\n") +
		  guides(fill=guide_legend(title=c("Copy number")))
		  
pdf(file=out_file[[6]], width=6, height=6)
print(plot.0)
dev.off()

plot.0 =  ggplot(mutation_summary, aes(x = VAF, y = DP, fill=total_cn)) +
		  geom_point(alpha=.85, size=2.5, shape=21) +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=8), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nVAF (%)\n", y="DP\n") +
		  scale_x_log10() +
		  annotation_logticks(side="b") +
		  coord_cartesian(xlim=c(5,100)) +
		  guides(fill=guide_legend(title=c("Copy number")))
		  
pdf(file=out_file[[7]], width=6, height=6)
print(plot.0)
dev.off()


tmp = mutation_summary %>%
	  group_by(cluster_id) %>%
	  summarize(
	  		n = n(),
	  		mean_ccf = mean(ccf),
	    	median_ccf = median(ccf),
	    	std_ccf = sd(ccf),
	    	min_ccf = min(ccf),
	    	max_ccf = max(ccf),
	    	mean_sd = mean(std),
	    	median_sd = median(std),
	    	std_sd = sd(std),
	    	min_sd = min(std),
	    	max_sd = max(std))
	    	
write_tsv(x=mutation_summary, path=out_file[[8]])
write_tsv(x=tmp, path=out_file[[9]])
