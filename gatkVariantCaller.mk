# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

##### DEFAULTS ######
REF ?= hg19
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))
LOGDIR = log/gatk.$(NOW)

##### OPTIONS ######
# HARD_FILTER_SNPS = true/false (default: true)
# 	Filter snps using hard thresholds instead of dynamically
# 	Dynamic thresholds require a lot of data (e.g. a genome-wide coverage)
# POOL_SNP_RECAL = true/false (default: false)
# 	Pool snps for score recalibration
# SPLIT_CHR = true/false (default: true)
# 	not advised for amplicon libraries which have low/no coverage in some chromosomes

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

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
VARIANTS_VCF := $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.vcf)
SNPS_VCFS := $(foreach sample,$(SAMPLES),vcf/$(sample).gatk_snps.eff.nsfp.vcf)
SNPS_TABLES := $(foreach sample,$(SAMPLES),tables/$(sample).gatk_snps.eff.nsfp.novel.txt)
INDELS_VCFS := $(foreach sample,$(SAMPLES),vcf/$(sample).gatk_indels.eff.nsfp.vcf)
INDELS_TABLES := $(foreach sample,$(SAMPLES),tables/$(sample).gatk_indels.eff.nsfp.novel.txt)

.PHONY : all variants_vcf variants_table recal_sample_bams realn_chr_bams realn_sample_bams variants_novel_table sample_variants clean reports all_variants
all : snps_vcf indels_vcf snps_table indels_table
variants_vcf : $(VARIANTS_VCF)
snps_vcf : $(SNPS_VCFS)
indels_vcf : $(INDELS_VCFS)
snps_table : $(SNPS_TABLES)
indels_table : $(INDELS_TABLES)
reports : gatk/reports/snps_filter.grp gatk/reports/indels_filter.grp gatk/reports/all.grp


include ~/share/modules/gatk.mk
