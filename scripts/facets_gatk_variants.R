#!/usr/bin/env Rscript

suppressMessages(pacman::p_load(dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,crayon))

#------
# input
#------

# read samples file
samples <-
	read.delim("sample_sets.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE) %>%
	filter(!grepl("#",V1)) %>%
	(function(sets){
		tumor <-
			apply(sets,1, function(row) row %>% list.filter (.!="") %>% head(-1)) %>% unlist
		normal <-
			apply(sets,1, function(row) row %>% list.filter (.!="") %>% tail(1)) %>%
			rep(apply(sets,1,function(row) row %>% list.filter(.!="") %>% length-1))
		data.frame(normal,tumor,stringsAsFactors=FALSE) %>%
		tbl_df
	})

#-------------------------
# create combined GATK vcf
#-------------------------

vcf.paths <-
	paste("ls ",samples$tumor %>% paste("gatk/vcf/",.,"_*.variants.snps.filtered.vcf",sep="",collapse=" ")) %>%
	system(intern=TRUE)

system("mkdir -p facets/gatk_variant_input &>/dev/null")

#filter vcfs for quality & depth
vcf.paths %>%
substr(10,(nchar(.)-4)) %>%
lapply(. %>% paste("vcftools --vcf gatk/vcf/",.,".vcf --minGQ 20 --minDP 8 --recode --out facets/gatk_variant_input/",.,sep="") %>% system)

# gzip files
list.files("facets/gatk_variant_input",pattern="vcf$") %>%
list.filter(.!="all.variants.snps.filtered.recode.vcf") %>%
paste("bgzip -c facets/gatk_variant_input/",.," > facets/gatk_variant_input/",.,".gz",sep="") %>%
lapply(. %>% system)

# tabix files
list.files("facets/gatk_variant_input",pattern="vcf.gz$") %>%
paste("tabix -p vcf facets/gatk_variant_input/",.,sep="") %>%
lapply(. %>% system)

# combine variants into single file & gzip
list.files("facets/gatk_variant_input",pattern="vcf.gz$") %>%
paste("facets/gatk_variant_input/",.,sep="",collapse=" ") %>%
paste("vcf-merge ",.," > facets/gatk_variant_input/all.variants.snps.filtered.recode.vcf",sep="") %>%
system

# gzip combined file
"bgzip -c facets/gatk_variant_input/all.variants.snps.filtered.recode.vcf > facets/gatk_variant_input/all.variants.snps.filtered.recode.vcf.gz" %>%
system