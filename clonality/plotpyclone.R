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

mutation_summary = read_tsv(file=paste0("sufam/", opt$sample_set, ".tsv")) %>%
				   mutate(mutation_id = paste0(Gene_Symbol, "_", HGVSp))
index = apply(mutation_summary[,paste0("DP_", tumor_samples)], 1, function(x) {sum(x>=min_depth)})==length(tumor_samples)
mutation_summary = mutation_summary[index,,drop=FALSE]
pyclone_summary = read_tsv(file=paste0("pyclone/", opt$sample_set, "/report/pyclone.tsv"), col_types = cols(.default = col_character())) %>%
				  type_convert() %>%
				  full_join(mutation_summary, by="mutation_id") %>%
				  arrange(cluster_id) %>%
				  mutate(mutation_type = ifelse(Variant_Caller=="mutect", "SNV", "Indel")) %>%
				  mutate(nref = nchar(Ref)) %>%
				  mutate(nalt = nchar(Alt)) %>%
				  filter(nref<=2 & nalt<=2)

df = pyclone_summary[,c("mutation_id", "cluster_id", "mutation_type"),drop=FALSE]
for (i in 1:length(tumor_samples)) {
	x = pyclone_summary[,tumor_samples[i]] %>%
		.[[1]]
	c_x = pyclone_summary %>%
		  .[[paste0("CALL_", tumor_samples[i])]]
	m_x = pyclone_summary %>%
		  .[[paste0("MAF_", tumor_samples[i])]]		  
	x[x<.025 | c_x==0 | m_x<.05] = 0
	df = cbind(df, x)
	colnames(df)[i+3] = tumor_samples[i]
}
index = apply(df[,tumor_samples], 1, function(x) {sum(x==0)})==length(tumor_samples)
df = df[!index,,drop=FALSE]
pyclone_summary = pyclone_summary[!index,,drop=FALSE]
index = apply(pyclone_summary[,paste0("DP_", tumor_samples)], 1, function(x) {sum(x>=500)})>=1
df = df[!index,,drop=FALSE]
pyclone_summary = pyclone_summary[!index,,drop=FALSE]


pyclone_summary[,tumor_samples] = df[,tumor_samples]


clusters = table(pyclone_summary$cluster_id)
if (any(clusters==1)) {
	pyclone_summary = pyclone_summary %>%
					  filter(!(cluster_id %in% names(clusters)[clusters==1]))
}

df = pyclone_summary[,c("mutation_id", "cluster_id", "mutation_type"),drop=FALSE]
for (i in 1:length(tumor_samples)) {
	x = pyclone_summary[,tumor_samples[i]] %>%
		.[[1]]
	df = cbind(df, x)
	colnames(df)[i+3] = tumor_samples[i]
}


pdf(file=paste0("pyclone/", opt$sample_set, "/report/pyclone.pdf"), width=6.5, height=6)
for (i in 1:(length(tumor_samples)-1)) {
	for (j in (i+1):length(tumor_samples)) {
		x = df[,tumor_samples[i]]
		y = df[,tumor_samples[j]]
		z1 = df[,"cluster_id"]
		z2 = df[,"mutation_type"]
		tmp.0 = data_frame(x=x, y=y, z1=factor(z1, ordered=TRUE), z2=z2)
		plot.0 =  ggplot(tmp.0, aes(x=x, y=y, fill=z1, color=z1, shape=z2)) +
				  geom_point(alpha = .55, size=2.5) +
				  theme_classic() +
				  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
				  theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
				  labs(x=paste0("\n",tumor_samples[i],"\n"), y=paste0("\n",tumor_samples[j],"\n")) +
				  guides(color=guide_legend(title=c("Cluster")), shape=guide_legend(title=c("Type"))) +
				  guides(fill=FALSE)
		print(plot.0)
	}
}
dev.off()

write_tsv(pyclone_summary, path=paste0("pyclone/", opt$sample_set, "/report/summary.tsv"))
