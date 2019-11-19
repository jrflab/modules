#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--input_file", default = NA, type = 'character', help = "file name and path"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_vars = read_tsv(file=opt$input_file, col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
all_tumors = all_vars %>%
			 .[["TUMOR_SAMPLE"]]

all_normals = all_vars %>%
			 .[["NORMAL_SAMPLE"]]

all_patients = unique(paste0(all_tumors, "_", all_normals))

all_vars = all_vars %>%
	   filter(Variant_Classification=="Frame_Shift_Del" | Variant_Classification=="In_Frame_Del") %>%
	   filter((grepl("varscan", variantCaller) & grepl("strelka", variantCaller)) |
	   	  ((grepl("platypus", variantCaller) & grepl("scalpel", variantCaller)) & ((nchar(REF)-nchar(ALT))>4) & Variant_Classification!="In_Frame_Del") |
	   	  ((grepl("platypus", variantCaller) & grepl("lancet", variantCaller)) & ((nchar(REF)-nchar(ALT))>4) & Variant_Classification!="In_Frame_Del"))

patient_summary = data_frame(SAMPLE_UUID = all_patients)
del_count = all_vars %>%
			mutate(SAMPLE_UUID = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
			dplyr::group_by(SAMPLE_UUID) %>%
			dplyr::summarize(del_count = n())
mean_delen = all_vars %>%
			 mutate(del_len = nchar(REF)) %>%
			 mutate(SAMPLE_UUID = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
			 dplyr::group_by(SAMPLE_UUID) %>%
			 dplyr::summarize(mean_delen = mean(del_len))
median_delen = all_vars %>%
			   mutate(del_len = nchar(REF)) %>%
			   mutate(SAMPLE_UUID = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
			   dplyr::group_by(SAMPLE_UUID) %>%
			   dplyr::summarize(median_delen = median(del_len))
deln4_count = all_vars %>%
			  mutate(del_len = nchar(REF)) %>%
			  mutate(SAMPLE_UUID = paste0(TUMOR_SAMPLE, "_", NORMAL_SAMPLE)) %>%
			  dplyr::group_by(SAMPLE_UUID) %>%
			  dplyr::summarize(deln4_count = sum(del_len>=4))
			  
'getSeqFrom' <- function(chr, start, end)
{
	ret = as.character(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, names=chr, start=start, end=end, strand="+", as.character=TRUE))
	return(invisible(ret))
}


'checkHomLen' <- function(deleted, next50)
{
	ret = 0
	for (i in 1:nchar(deleted)) {
		if (substr(deleted, 1, i) == substr(next50, 1, i)) {
			ret = i
		}
	}
    return(invisible(ret))
}

hml_down = hml_up = NULL
for (i in 1:nrow(all_vars)) {
	chr = paste0("chr", all_vars[i,"CHROM"])
	start = as.numeric(all_vars[i,"POS"])
	n = as.numeric(nchar(all_vars[i,"REF"]))-1
	
	deleted = getSeqFrom(chr = chr, start = start, end = start + n)
	prevn = getSeqFrom(chr = chr, start = start - n - 1, end = start - 1)
	nextn = getSeqFrom(chr = chr, start = start + n + 1, end = start + 2*n + 1)
	
	hml_down = c(hml_down, checkHomLen(deleted = deleted, next50 = prevn))
	hml_up = c(hml_up, checkHomLen(deleted = deleted, next50 = nextn))
}

mh_3 = data_frame(SAMPLE_UUID = paste0(all_vars$TUMOR_SAMPLE, "_", all_vars$NORMAL_SAMPLE),
				   del_len = nchar(all_vars$REF),
				   max_mhlen_5p = hml_down,
				   max_mhlen_3p = hml_up,
	  			   max_mhlen = apply(cbind(hml_down, hml_up), 1, max)) %>%
	    filter(del_len >= 4) %>%
	    mutate(is_3 = ifelse(max_mhlen>=3, 1, 0)) %>%
	    dplyr::group_by(SAMPLE_UUID) %>%
	    dplyr::summarize(deln4_mhlen_3_counts = sum(is_3))

mhl_3 = data_frame(SAMPLE_UUID = paste0(all_vars$TUMOR_SAMPLE, "_", all_vars$NORMAL_SAMPLE),
				   del_len = nchar(all_vars$REF),
				   max_mhlen_5p = hml_down,
				   max_mhlen_3p = hml_up,
	  			   max_mhlen = apply(cbind(hml_down, hml_up), 1, max)) %>%
	    filter(del_len >= 4) %>%
	    filter(max_mhlen >= 3) %>%
	    dplyr::group_by(SAMPLE_UUID) %>%
	    dplyr::summarize(deln4_mhlen_3_avg_deln = mean(del_len))				 

patient_summary = left_join(patient_summary, del_count, by="SAMPLE_UUID") %>%
				  left_join(mean_delen, by="SAMPLE_UUID") %>%
				  left_join(median_delen, by="SAMPLE_UUID") %>%
				  left_join(deln4_count, by="SAMPLE_UUID") %>%
				  left_join(mh_3, by="SAMPLE_UUID") %>%
				  left_join(mhl_3, by="SAMPLE_UUID") %>%
				  mutate(delmh_prop = deln4_mhlen_3_counts/del_count) %>%
				  mutate(delmh_del4n_prop = deln4_mhlen_3_counts/deln4_count)
				  
write_tsv(patient_summary, path="summary/tsv/delmh_summary.tsv")
