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

in_file = paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/pyclone.tsv")
out_file = list(
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_std.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/report/histogram_ccf.pdf")
)

mutation_summary = read_tsv(file=in_file, col_types = cols(.default = col_character()), col_names = c("ID", "CCF", "STD", "C_ID")) %>%
				   type_convert() %>%
				   mutate(C_ID = factor(C_ID)) %>%
				   mutate(CCF = as.numeric(CCF)) %>%
				   mutate(STD = as.numeric(STD)) %>%
				   slice(-1)

plot.0 =  ggplot(mutation_summary, aes(x=STD, color=C_ID)) +
		  geom_histogram(fill="salmon", alpha=0.5, position="identity") +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nSigma\n", y="Frequency\n") +
		  guides(color=guide_legend(title=c("Cluster")))
pdf(file=out_file[[1]], width=6, height=6)
print(plot.0)
dev.off()
		 
plot.0 =  ggplot(mutation_summary, aes(x=CCF, color=C_ID)) +
		  geom_histogram(fill="salmon", alpha=0.5, position="identity") +
		  theme_classic() +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nCCF\n", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  guides(color=guide_legend(title=c("Cluster")))
pdf(file=out_file[[2]], width=6, height=6)
print(plot.0)
dev.off()
