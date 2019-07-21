#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(
					make_option("--sample_name", default = NA, type = 'character', help = "tumor sample name")
				  )
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

mutation_summary = read_tsv(file="summary/tsv/mutation_summary.tsv", col_types = cols(.default = col_character()))  %>%
 				   type_convert() %>%
 				   filter(variantCaller=="mutect") %>%
 				   filter(TUMOR_SAMPLE==opt$sample_name) %>%
 				   mutate(CHROM = paste0("chr", CHROM)) %>%
 				   select(sample_id = TUMOR_SAMPLE, chrom=CHROM, pos=POS, ref=REF, alt=ALT)

signature_input = mut.to.sigs.input(mut.ref = data.frame(mutation_summary),
									sample.id = "sample_id", 
									chr = "chrom", 
									pos = "pos", 
									ref = "ref", 
									alt = "alt")
									
extracted_signatures = whichSignatures(tumor.ref = signature_input,
									   signatures.ref = signatures.cosmic,
									   contexts.needed = TRUE)
									   
save(list=ls(all=TRUE), file=paste0("deconstructsigs/signatures/", opt$sample_name, ".RData"))
