#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_name", default = NA, type = 'character', help = "tumor sample name")
				  )
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

load(file=paste0("deconstructsigs/signatures/", opt$sample_name, ".RData"))

df = data_frame(percentage = 100*as.vector(extracted_signatures$tumor),
				trinucleotide_context = colnames(extracted_signatures$tumor)) %>%
				mutate(ref = rep(c("C", "T"), each=48)) %>%
				mutate(alt = rep(c("A", "G", "T", "A", "C", "G"), each=16)) %>%
				mutate(base_change = factor(paste0(ref, ">", alt)))
				

plot.0 = ggplot(df, aes(x=trinucleotide_context, y=percentage, fill=base_change)) +
		 geom_bar(stat="identity") +
		 facet_wrap(~base_change, ncol = 6, nrow = 1, scales = "free_x") +
  		 ylab("\nFraction (%)\n") +
  		 xlab(" ") +
  		 theme_bw(base_size=15) +
  		 theme(axis.text.y = element_text(size=14), axis.text.x = element_text(size=10, angle=90), legend.position="none")


pdf(file=paste0("deconstructsigs/plots/trint_context/", opt$sample_name, ".pdf"), width=25, height=5)
print(plot.0)
dev.off()

palette = colorRampPalette(brewer.pal(9, "Set1"))
cols = palette(30)
names(cols) = 1:30

df = data.frame(percentage = 100*as.numeric(extracted_signatures$weights[1,]),
				signature_name = colnames(extracted_signatures$weights)) %>%
				mutate(signature_name = gsub(pattern="Signature.", replacement="", signature_name)) %>%
				filter(percentage!=0) %>%
				arrange(signature_name) %>%
				mutate(signature_name = factor(signature_name)) %>%
				mutate(lab.ypos = cumsum(percentage) - 0.5*percentage)
				
plot.0  = ggplot(df, aes(x = "", y = percentage, fill = signature_name)) +
		  geom_bar(width = 1, stat = "identity", color = "white") +
		  scale_fill_manual(values=cols) +
		  coord_polar("y", start = 0) +
		  geom_text(aes(y = lab.ypos, label = signif(percentage,3)), color = "white") +
		  guides(fill=guide_legend(title="Signature")) +
		  theme_void()
		  
pdf(file=paste0("deconstructsigs/plots/signature_exposures/", opt$sample_name, ".pdf"), width=6, height=6)
print(plot.0)
dev.off()
