#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("signature.tools.lib"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--option", default = NA, type = 'character', help = "type of analysis"),
		  make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$option) == 1) {
	vcf = readr::read_tsv(file = "summary/tsv/all.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      dplyr::filter(CHROM %in% c(1:22, "X")) %>%
	      dplyr::mutate(CHROM = case_when(
		      CHROM == "X" ~ "23",
		      TRUE ~ CHROM
	      )) %>%
	      readr::type_convert() %>%
	      dplyr::arrange(CHROM, POS) %>%
	      dplyr::mutate(TUMOR_NORMAL = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
	      dplyr::filter(TUMOR_NORMAL == as.character(opt$sample_name)) %>%
	      dplyr::filter(variantCaller == "mutect") %>%
	      dplyr::filter(TUMOR_DP>=10 & NORMAL_DP>=10) %>%
	      dplyr::mutate(CHROM = as.character(CHROM)) %>%
	      dplyr::mutate(CHROM = ifelse(CHROM == "23", "X", CHROM)) %>%
	      dplyr::select(`#CHROM` = CHROM,
			    POS = POS,
			    ID = ID,
			    REF = REF,
			    ALT = ALT)
	cat("##fileformat=VCFv4.1\n", file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), append = FALSE)
	readr::write_tsv(x = vcf, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".snv.vcf"), col_names = TRUE, append = TRUE)

} else if (as.numeric(opt$option) == 2) {
	vcf = readr::read_tsv(file = "summary/tsv/all.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      dplyr::filter(CHROM %in% c(1:22, "X")) %>%
	      dplyr::mutate(CHROM = case_when(
		      CHROM == "X" ~ "23",
		      TRUE ~ CHROM
	      )) %>%
	      readr::type_convert() %>%
	      dplyr::arrange(CHROM, POS) %>%
	      readr::type_convert() %>%
	      dplyr::mutate(TUMOR_NORMAL = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
	      dplyr::filter(TUMOR_NORMAL == as.character(opt$sample_name)) %>%
	      dplyr::filter(grepl("varscan", variantCaller, fixed = TRUE)) %>%
	      dplyr::filter(grepl("strelka", variantCaller, fixed = TRUE)) %>%
	      dplyr::filter(TUMOR_DP>=10 & NORMAL_DP>=10) %>%
	      dplyr::mutate(CHROM = as.character(CHROM)) %>%
	      dplyr::mutate(CHROM = ifelse(CHROM == "23", "X", CHROM)) %>%
	      dplyr::select(`#CHROM` = CHROM,
			    POS = POS,
			    ID = ID,
			    REF = REF,
			    ALT = ALT)
	cat("##fileformat=VCFv4.1\n", file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".indel.vcf"), append = FALSE)
	readr::write_tsv(x = vcf, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".indel.vcf"), col_names = TRUE, append = TRUE)

} else if (as.numeric(opt$option) == 3) {
	cn = readr::read_tsv(file = paste0("facets/cncf/", as.character(opt$sample_name), ".txt"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
     	     dplyr::mutate(chrom = as.character(chrom)) %>%
	     dplyr::mutate(chrom = ifelse(chrom == "23", "X", chrom)) %>%
	     dplyr::mutate(seg_no = seg,
			   Chromosome = chrom,
			   chromStart = loc.start,
			   chromEnd = loc.end,
			   total.copy.number.inNormal = 2,
			   minor.copy.number.inNormal = 1,
			   total.copy.number.inTumour = tcn.em,
			   minor.copy.number.inTumour = lcn.em) %>%
	     dplyr::mutate(total.copy.number.inTumour = case_when(
		     		is.na(total.copy.number.inTumour) ~ 2,
		     		TRUE ~ total.copy.number.inTumour
	     )) %>%
	     dplyr::mutate(minor.copy.number.inTumour = case_when(
		     		is.na(minor.copy.number.inTumour) ~ 2,
		     		TRUE ~ minor.copy.number.inTumour
	     )) %>%
	     dplyr::select(seg_no,
			   Chromosome,
			   chromStart,
			   chromEnd,
			   total.copy.number.inNormal,
			   minor.copy.number.inNormal,
			   total.copy.number.inTumour,
			   minor.copy.number.inTumour)
	     
	readr::write_tsv(x = cn, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".cn.txt"), col_names = TRUE, append = FALSE)
	
} else if (as.numeric(opt$option) == 4) {
	sv = readr::read_tsv(file = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".merged.bedpe"), col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
	     dplyr::mutate(svclass = case_when(
		     svclass == "BND" ~ "translocation",
		     svclass == "DEL" ~ "deletion",
		     svclass == "DUP" ~ "tandem-duplication",
		     svclass == "INS" ~ "insertion",
		     svclass == "INV" ~ "inversion",
	     	     TRUE ~ svclass
	     )) %>%
             dplyr::filter(svclass != "inversion") %>%
     	     dplyr::mutate(sample = as.character(opt$sample_name))
	     
	readr::write_tsv(x = sv, path = paste0("hr_detect/", as.character(opt$sample_name), "/", as.character(opt$sample_name), ".sv.bedpe"), col_names = TRUE, append = FALSE)

	
} else if (as.numeric(opt$option) == 5) {
	sample_names = unlist(strsplit(x = as.character(opt$sample_name), split = " ", fixed = TRUE))
	snv_files = unlist(lapply(sample_names, function(x) { paste0("hr_detect/", x, "/", x, ".snv.vcf") }))
	indel_files = unlist(lapply(sample_names, function(x) { paste0("hr_detect/", x, "/", x, ".indel.vcf.bgz") }))
	cn_files = unlist(lapply(sample_names, function(x) { paste0("hr_detect/", x, "/", x, ".cn.txt") }))
	sv_files = unlist(lapply(sample_names, function(x) { paste0("hr_detect/", x, "/", x, ".sv.bedpe") }))
	
	names(snv_files) = names(indel_files) = names(cn_files) = names(sv_files) <- sample_names
	
	snv_cat_list = list()
	for (i in 1:length(snv_files)) {
		snv_tab = readr::read_tsv(file = snv_files[i], col_names = TRUE, comment = "##", col_types = cols(.default = col_character())) %>%
			  readr::type_convert() %>%
			  dplyr::select(chr = `#CHROM`,
				        position = POS,
				        REF,
					ALT) %>%
			 as.data.frame()
		res = tabToSNVcatalogue(subs = snv_tab, genome.v = "hg19")
		colnames(res$catalogue) = sample_names[i]
		snv_cat_list[[i]] = res$catalogue
	}
	snv_catalogues = do.call(cbind, snv_cat_list)
	
	sigsToUse = c(1,2,3,5,6,8,13,17,18,20,26,30)
	subs_fit_res = Fit(catalogues = snv_catalogues,
			   signatures = COSMIC30_subs_signatures[,sigsToUse],
			   useBootstrap = TRUE,
			   nboot = 100,
			   nparallel = 4)
	snv_exp = subs_fit_res$exposures
	
	col_hrdetect = c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
	input_matrix = matrix(NA, nrow = length(sample_names), ncol = length(col_hrdetect), dimnames = list(sample_names, col_hrdetect))
	input_matrix[rownames(snv_exp),"SNV3"] = snv_exp[,"Signature3"]
	input_matrix[rownames(snv_exp),"SNV8"] = snv_exp[,"Signature8"]
	res =  HRDetect_pipeline(input_matrix,
				 genome.v = "hg19",
				 SNV_signature_version = "COSMICv2",
				 SV_bedpe_files = sv_files,
				 Indels_tab_files = indel_files,
				 CNV_tab_files = cn_files,
				 nparallel = 4)
	
	readr::write_tsv(x = res$hrdetect_output %>%
			     dplyr::as_tibble() %>%
			     dplyr::mutate(sample_name = sample_names),
			path = "hr_detect/summary.txt", append = FALSE, col_names = TRUE)
}
