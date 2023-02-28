#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))

options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(make_option("--input_file", default = 'snp_vcf/snps_ft.vcf', help = "input file"),
		make_option("--output_file", default = 'snp_vcf/snps_ft.pdf', help = "output file"),
		make_option("--sample_pairs", default = NA, help = "sample pairs"),
		make_option("--genome", default = 'b37', help = "genome build"))

parser <- OptionParser(usage = "%prog vcf.files", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

vcf = readVcf(as.character(opt$input_file), as.character(opt$genome))
gt = geno(vcf)$GT
ad = geno(vcf)$AD
af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
X[is.na(af)] = NA
X[af > 0.15 & af < 0.95] = 1
X[af >= 0.95] = 2
X[!gt %in% c("0/0", "0/1", "1/1")] = NA

gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
dt = as.matrix(dist(t(gt)))
		      
tumor_samples = unlist(lapply(strsplit(x = unlist(strsplit(x = as.character(opt$sample_pairs), split = " ")), split = "_"), function(x) { x[1] }))
normal_samples = unlist(lapply(strsplit(x = unlist(strsplit(x = as.character(opt$sample_pairs), split = " ")), split = "_"), function(x) { x[2] }))
sample_pairs = dplyr::tibble(tumor_samples = factor(c(tumor_samples, unique(normal_samples)), levels = rownames(dt), ordered = TRUE),
			     normal_samples = c(normal_samples, unique(normal_samples))) %>%
	       dplyr::arrange(tumor_samples) %>%
	       dplyr::mutate(normal_samples = factor(normal_samples, levels = unique(normal_samples), ordered = TRUE))
cluster_color = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(sample_pairs %>% .[["normal_samples"]])))
names(cluster_color) = sort(unique(sample_pairs %>% .[["normal_samples"]]))
row_annot = rowAnnotation(
	cluster_id = sample_pairs %>% .[["normal_samples"]],
	col = list(cluster_id = cluster_color),
		   show_annotation_name = FALSE,
		   simple_anno_size = unit(.5, "cm"),
		   show_legend = FALSE
)
col_annot = columnAnnotation(
	cluster_id = sample_pairs %>% .[["normal_samples"]],
	col = list(cluster_id = cluster_color),
		   show_annotation_name = FALSE,
		   simple_anno_size = unit(.5, "cm"),
		   show_legend = FALSE
)
col_pal = c(rep("#662506", 3),
	    rev(brewer.pal(n = 7, name = "YlOrBr")),
	    rep("#fff7bc", 3))

pdf(as.character(opt$output_file), height = 21, width = 22)
draw(Heatmap(matrix = dt,
	     name = " ",
	     rect_gp = gpar(col = "white"),
	     border = NA,
	     col = col_pal,
	     cluster_rows = TRUE,
	     show_row_dend = TRUE,
	     row_dend_width = unit(3, "cm"),
	     row_names_side = "right",
	     row_names_gp = gpar(fontsize = 12),
	     show_row_names = TRUE,
	     left_annotation = row_annot,
	          
	     show_column_names = TRUE,
	     column_names_side = "bottom",
	     column_names_gp = gpar(fontsize = 12),
	     cluster_columns = TRUE,
	     show_column_dend = TRUE,
	     column_dend_height = unit(3, "cm"),
	     top_annotation = col_annot,
	     
	     use_raster = FALSE,
	     show_heatmap_legend = TRUE,
	     heatmap_legend_param = list(legend_height = unit(5, "cm"), legend_width = unit(5, "cm"))))
dev.off()
