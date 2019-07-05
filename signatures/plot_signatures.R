#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("Palimpsest"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_name", default = NA, type = 'character', help = "tumor sample name")
				  )
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

'plot96_mutation_spectrum' <- function (vcf, sample.col = "sample", mutcat3.col = "mutcat3",
										ymax = NULL, averageProp = FALSE, plot.file = NULL)
{
    bases <- c("A", "C", "G", "T")
    ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
    mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
    types96 <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
    types96 <- sapply(types96, function(z) {
        sub("\\.", substr(z, 1, 1), z)
    })
    context <- substr(types96, 4, 6)
    nsamp <- length(unique(vcf[, sample.col]))
    if (averageProp & nsamp > 1) {
        tmp <- makeMutypeMatFromVcf(vcf, sample.col = "CHCID", 
            mutcat.col = "mutcat3", mutypes = types96)
        freq <- apply(tmp, 1, mean)
    }
    else {
        freq <- sapply(types96, function(z) {
            mean(vcf[, mutcat3.col] == z, na.rm = T)
        })
    }
    if (!is.null(plot.file)) {
        pdf(plot.file, width = 24, height = 5)
    }
    col96 <- c(rep("skyblue3", 16), rep("black", 16), rep("red", 
        16), rep("grey", 16), rep("green", 16), rep("pink", 16))
    labs <- c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16), 
        rep("T>A", 16), rep("T>C", 16), rep("T>G", 16))
    if (is.null(ymax)) {
        ymax <- 100*ceiling(max(freq) * 100)/100
        ymax <- ifelse(ymax>10, 30, 10)
    }
    bp <- barplot(freq*100, col = col96, border = col96, las = 2, 
        width = 1, space = .35, yaxt = "n", xaxt = "n", ylim = c(0, 
            ymax * 1.2))
    title(ylab = "Fraction of mutations (%)", mgp = c(1, 1, 0), 
        cex.lab = 1.6)
    axis(1, at = bp, labels = context, pos = 0, las = 2, cex.axis = 1.5, 
        tick = F, cex.axis = 1, lwd=-1)
    if (ymax==40) {
	    axis(2, at = c(0,10,20,30,40), labels=c(0,10,20,30,40), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==30) {
	    axis(2, at = c(0,5,10,15,20,25,30), labels=c(0,5,10,15,20,25,30), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==20) {
		axis(2, at = c(0,5,10,15,20), labels=c(0,5,10,15,20), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==10) {
		axis(2, at = c(0,2,4,6,8,10), labels=c(0,2,4,6,8,10), pos = 0, las = 1, cex.axis = 1.5)
	}
    for (i in seq(1, 81, by = 16)) {
        rect(bp[i], par()$usr[4], bp[i + 15], par()$usr[4] - 
            0.05 * diff(par()$usr[3:4]), col = col96[i], border = col96[i])
        text((bp[i] + bp[i + 15])/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), 
            labels = labs[i], xpd = TRUE, cex = 2)
    }
    if (!is.null(plot.file)) {
        dev.off()
    }
}

load(file=paste0("deconstructsigs/signatures/", opt$sample_name, ".RData"))

## barplot of base changes with 3' and 5' context
colnames(mutation_summary) = c("Sample", "CHROM", "POS", "REF", "ALT")
mutation_summary = cbind(mutation_summary, "Type"=rep("SNV", nrow(mutation_summary)))
vcf = preprocessInput_snv(input_data = mutation_summary,
                          ensgene = ensgene,
                          reference_genome = BSgenome.Hsapiens.UCSC.hg19)
patient_ids = unique(vcf$Sample)
pdf(file=paste0("deconstructsigs/plots/context/", opt$sample_name, ".pdf"), width=18, height=5)
plot96_mutation_spectrum(vcf, ymax=20, sample.col = "Sample",  plot.file = NULL)
dev.off()

## pie-charts of signatures
palette = colorRampPalette(brewer.pal(9, "Set1"))
cols = palette(30)
names(cols) = 1:30

df = data.frame(percentage = 100*as.numeric(extracted_signatures$weights[1,]),
				signature_name = colnames(extracted_signatures$weights)) %>%
				mutate(signature_name = as.numeric(gsub(pattern="Signature.", replacement="", signature_name))) %>%
				arrange(signature_name) %>%
				filter(percentage!=0) %>%
				mutate(signature_name = factor(signature_name, ordered=TRUE, levels=sort(signature_name))) %>%
				mutate(lab.ypos = cumsum(percentage) - 0.5*percentage)
				
plot.0  = ggplot(df, aes(x = "", y = percentage, fill = signature_name)) +
		  geom_bar(width = 1, stat = "identity", color = "white") +
		  scale_fill_manual(values=cols) +
		  coord_polar("y", start = 0) +
		  geom_text(aes(y = lab.ypos, label = paste0(signif(percentage,3), "%")), color = "white") +
		  guides(fill=guide_legend(title="Signature")) +
		  theme_void()
		  
pdf(file=paste0("deconstructsigs/plots/exposures/", opt$sample_name, ".pdf"), width=6, height=6)
print(plot.0)
dev.off()
