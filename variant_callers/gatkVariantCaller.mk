# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

##### DEFAULTS ######
LOGDIR ?= log/gatk.$(NOW)

##### OPTIONS ######
# HARD_FILTER_SNPS = true/false (default: true)
# 	Filter snps using hard thresholds instead of dynamically
# 	Dynamic thresholds require a lot of data (e.g. a genome-wide coverage)
# POOL_SNP_RECAL = true/false (default: false)
# 	Pool snps for score recalibration
# SPLIT_CHR = true/false (default: true)
# 	not advised for amplicon libraries which have low/no coverage in some chromosomes

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

VPATH ?= bam

.DELETE_ON_ERROR:

.SECONDARY: 

######## BUILD ORDER ##########
# determine realignment targets (per chromosome if SPLIT_CHR=true (default))
# realignment of bams (per chromosome if SPLIT_CHR=true)
# recalibration of base qualities (per chromosome if SPLIT_CHR=true)
# reduction of bam size (per chromosome if SPLIT_CHR=true)
# if (SPLIT_CHR=true) merge realigned, recalibrated, reduced chromosomes bams
# predict variants 
# split indels and snps
# annotate using snpEff
# if (HARD_FILTER_SNPS = true) hard snp filterings
# 	else if (POOL_SNP_RECAL = true) pool calls across samples for recalibration
# 		else calibrated snp filtering per sample
# hard filtering of indels
# create indel and snp tables
# create novel indel/snp tables

##### MAIN TARGETS ######
VARIANT_TYPES = gatk_snps gatk_indels

PHONY += gatk gatk_vcfs gatk_mafs

gatk : gatk_vcfs #gatk_mafs reports

gatk_vcfs : $(foreach type,$(VARIANT_TYPES),$(foreach sample,$(SAMPLES),vcf/$(sample).$(type).vcf))
gatk_mafs : $(foreach type,$(VARIANT_TYPES),$(foreach sample,$(SAMPLES),maf/$(sample).$(type).maf))


include modules/variant_callers/gatk.mk

.PHONY : $(PHONY)
