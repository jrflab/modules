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

file_paths = list(
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/pyclone.tsv"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/hist_std.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/hist_ccf.pdf"),
	paste0("pyclone/", tumor_sample, "_", normal_sample, "/scatter_ccf_std.pdf")
)

mutation_summary = read_tsv(file=file_paths[[1]], col_types = cols(.default = col_character()), col_names = c("ID", "CCF", "STD", "C_ID")) %>%
				   type_convert() %>%
				   mutate(FACET_1 = "Standard Deviation") %>%
				   mutate(FACET_2 = "Cancer Cell Fraction") %>%
				   mutate(FACET_3 = "Std versus CCF") %>%
				   mutate(C_ID = factor(C_ID)) %>%
				   mutate(CCF = as.numeric(CCF)) %>%
				   mutate(STD = as.numeric(STD)) %>%
				   slice(-1)

plot.0 =  ggplot(mutation_summary, aes(x=STD, color=C_ID)) +
		  geom_histogram(fill="white", alpha=0.5, position="identity") +
		  theme_bw(base_size=15) +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nSigma\n", y="Frequency\n") +
		  theme(legend.justification = c(1, 0),
		 	    legend.position = c(1, .5),
		 	    legend.title = element_blank(),
		 	    legend.background = element_blank(),
		 	    legend.text=element_text(size=8)) +
		 facet_wrap(~FACET_1)
pdf(file=file_paths[[2]], width=6, height=6)
print(plot.0)
dev.off()
		 
		 
plot.0 =  ggplot(mutation_summary, aes(x=CCF, color=C_ID)) +
		  geom_histogram(fill="white", alpha=0.5, position="identity") +
		  theme_bw(base_size=15) +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nCCF\n", y="Frequency\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  theme(legend.justification = c(1, 0),
		 	    legend.position = c(1, .5),
		 	    legend.title = element_blank(),
		 	    legend.background = element_blank(),
		 	    legend.text=element_text(size=8)) +
		 facet_wrap(~FACET_2)
pdf(file=file_paths[[3]], width=6, height=6)
print(plot.0)
dev.off()
				   
plot.0 =  ggplot(mutation_summary, aes(x=CCF, y=STD, fill=C_ID)) +
		  geom_point(alpha=.85, size=2.5, shape=21) +
		  theme_bw(base_size=15) +
		  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		  labs(x="\nCCF\n", y="Sigma\n") +
		  coord_cartesian(xlim=c(0,1)) +
		  theme(legend.justification = c(1, 0),
		 	    legend.position = c(1, .5),
		 	    legend.title = element_blank(),
		 	    legend.background = element_blank(),
		 	    legend.text=element_text(size=8)) +
		 facet_wrap(~FACET_3)
pdf(file=file_paths[[4]], width=6, height=6)
print(plot.0)
dev.off()



