#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("copynumber"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "type of analysis"),
		  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

'plot_log2_ratio' <- function(x)
{
   	par(mar=c(5, 5, 4, 2)+.1)
	plot(x = x$position, y = x$log2, type = "p", pch = ".", cex = 1, col = "grey75", axes = FALSE, frame = FALSE, xlab = "", ylab = "", main = "", ylim = c(-4,5))
	y = x %>%
	    dplyr::group_by(chromosome) %>%
	    dplyr::summarize(start = min(start_chr),
			     end = max(end_chr)) %>%
	    dplyr::mutate(chromosome = factor(chromosome, levels = c(1:22, "X"), ordered = TRUE)) %>%
	    dplyr::arrange(chromosome)
	points(x = c(y$start[1]-1E9, y$end[nrow(y)]), y = c(0, 0), type = "l", col = "grey20", lwd = 1.15)
	axis(1, at = c(y$start, y$end[nrow(y)]), labels = rep(" ", nrow(y)+1), cex.axis = 0.85, las = 1, tck = .035)
	axis(1, at = .5*(y$start + y$end), labels = y$chromosome, cex.axis = 0.85, las = 1)
	axis(2, at = c(-2, -1, 0, 1, 2), labels = c(-2, -1, 0, 1, 2), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio   "), line = 3.15, cex = 1)
}

'add_segmented' <- function(x)
{
	for (i in 1:nrow(x)) {
		points(x = c(x$Start_Position[i], x$End_Position[i]), y = rep(x$Log2_Ratio[i], 2), type = "l", col = "#e41a1c", lwd = 2.75)
	}
}

if (as.numeric(opt$option) == 1) {
	data = readr::read_tsv(file = paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::filter(weight>.1) %>%
	       dplyr::filter(chromosome != "Y") %>%
	       dplyr::mutate(chromosome = factor(chromosome, levels = c(1:22, "X"), ordered = TRUE))
	cytoband = data %>%
		   dplyr::group_by(chromosome) %>%
		   dplyr::summarize(start = min(start),
				    end = max(end)) %>%
		   dplyr::mutate(chromosome = factor(chromosome, levels = c(1:22, "X"), ordered = TRUE)) %>%
		   dplyr::arrange(chromosome) %>%
		   dplyr::mutate(end = cumsum(end))
	start = rep(0, nrow(cytoband))
	start[2:nrow(cytoband)] = cytoband$end[1:(nrow(cytoband)-1)] + cytoband$start[2:nrow(cytoband)]
	cytoband$start = start
	data = data %>%
	       dplyr::left_join(cytoband %>%
				dplyr::rename(start_chr = start,
					      end_chr = end),
			        by = "chromosome") %>%
	       dplyr::mutate(start = start + start_chr,
			     end = end + start_chr) %>%
	       dplyr::mutate(position = .5*(start + end)) %>%
	       dplyr::mutate(log2 = case_when(
		       log2 > 6 ~ 0,
		       log2 < (-4) ~ 0,
		       TRUE ~ log2
	       ))
	
	pdf(file = paste0("cnvkit/plots/log2/", opt$sample_name, ".pdf"), width = 8, height = 3.75)
	plot_log2_ratio(x = data)
	dev.off()

} else if (as.numeric(opt$option) == 2) {
	data = readr::read_tsv(file = paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::filter(weight>.1) %>%
	       dplyr::filter(chromosome != "Y")
	smoothed = winsorize(data = data %>% dplyr::select(chromosome, start, log2) %>% data.frame(), method = "mad")
	segmented = pcf(data = smoothed, kmin = 25, gamma = 50, normalize = TRUE, fast = FALSE) %>%
		    dplyr::as_tibble() %>%
		    dplyr::select(Sample_Name = sampleID, Chromosome = chrom, Arm = arm,
				  Start_Position = start.pos, End_Position = end.pos,
				  N = n.probes, Log2_Ratio = mean) %>%
		    dplyr::mutate(Sample_Name = opt$sample_name)
	readr::write_tsv(x = segmented, file = paste0("cnvkit/segmented/", opt$sample_name, ".txt"), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option) == 3) {
	data = readr::read_tsv(file = paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	       readr::type_convert() %>%
	       dplyr::filter(weight>.1) %>%
	       dplyr::filter(chromosome != "Y")
	smoothed = winsorize(data = data %>% dplyr::select(chromosome, start, log2) %>% data.frame(), method = "mad")
	segmented = pcf(data = smoothed, kmin = 25, gamma = 50, normalize = TRUE, fast = FALSE) %>%
		    dplyr::as_tibble() %>%
		    dplyr::select(Sample_Name = sampleID, Chromosome = chrom, Arm = arm,
				  Start_Position = start.pos, End_Position = end.pos,
				  N = n.probes, Log2_Ratio = mean) %>%
		    dplyr::mutate(Sample_Name = opt$sample_name)
	cytoband = data %>%
		   dplyr::group_by(chromosome) %>%
		   dplyr::summarize(start = min(start),
				    end = max(end)) %>%
		   dplyr::mutate(chromosome = factor(chromosome, levels = c(1:22, "X"), ordered = TRUE)) %>%
		   dplyr::arrange(chromosome) %>%
		   dplyr::mutate(end = cumsum(end))
	start = rep(0, nrow(cytoband))
	start[2:nrow(cytoband)] = cytoband$end[1:(nrow(cytoband)-1)] + cytoband$start[2:nrow(cytoband)]
	cytoband$start = start
	data = data %>%
	       dplyr::left_join(cytoband %>%
				dplyr::rename(start_chr = start,
					      end_chr = end),
			        by = "chromosome") %>%
	       dplyr::mutate(start = start + start_chr,
			     end = end + start_chr) %>%
	       dplyr::mutate(position = start) %>%
	       dplyr::mutate(log2 = case_when(
		       log2 > 6 ~ 0,
		       log2 < (-4) ~ 0,
		       TRUE ~ log2
	       ))
	segmented = segmented %>%
		    dplyr::left_join(cytoband %>%
				     dplyr::rename(Chromosome = chromosome,
						   start_chr = start,
					      	   end_chr = end),
			             by = "Chromosome") %>%
	       	    dplyr::mutate(Start_Position = Start_Position + start_chr,
				  End_Position = End_Position + start_chr)
	
	pdf(file = paste0("cnvkit/plots/segmented/", opt$sample_name, ".pdf"), width = 8, height = 3.75)
	plot_log2_ratio(x = data)
	add_segmented(x = segmented)
	dev.off()
}
