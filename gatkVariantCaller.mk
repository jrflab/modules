# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

##### DEFAULTS ######
REF ?= hg19
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))
LOGDIR = log/gatk.$(NOW)

VCF_GEN_IDS = GT AD GQ PL

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
ANN_TYPES = eff # annotated
EFF_TYPES = silent missense nonsilent_cds nonsilent
VARIANT_TYPES = gatk_snps gatk_indels

FILTER_SUFFIX := dp_ft.dbsnp
ifdef NORMAL_VCF
FILTER_SUFFIX := nft.$(FILTER_SUFFIX)
endif
FILTER_SUFFIX.gatk_snps := $(FILTER_SUFFIX).nsfp.chasm.fathmm
FILTER_SUFFIX.gatk_indels := $(FILTER_SUFFIX)
VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(type).$(FILTER_SUFFIX.$(type)).$(ann)))
TABLE_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(foreach eff,$(EFF_TYPES),$(type).$(FILTER_SUFFIX.$(type)).$(ann).tab.$(eff).pass.novel)))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(VCF_SUFFIXES),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(TABLE_SUFFIXES),tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(TABLE_SUFFIXES),tables/all.$(suff).txt)

ifdef SAMPLE_SETS
SS_FILTER_SUFFIX := dp_ft.som_ft.dbsnp
SS_FILTER_SUFFIX.gatk_snps := $(SS_FILTER_SUFFIX).nsfp.chasm.fathmm
SS_FILTER_SUFFIX.gatk_indels := $(SS_FILTER_SUFFIX)
SS_VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(type).$(SS_FILTER_SUFFIX.$(type)).$(ann)))
SS_TABLE_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(foreach ann,$(ANN_TYPES),$(foreach eff,$(EFF_TYPES),$(type).$(SS_FILTER_SUFFIX.$(type)).$(ann).tab.$(eff).pass.novel)))
VCFS += $(foreach set,$(SAMPLE_SETS),$(foreach suff,$(SS_VCF_SUFFIXES),vcf/$(set).$(suff).vcf))
TABLES += $(foreach set,$(SAMPLE_SETS),$(foreach suff,$(SS_TABLE_SUFFIXES),tables/$(set).$(suff).txt))
TABLES += $(foreach suff,$(SS_TABLE_SUFFIXES),tables/allSS.$(suff).txt)
endif

.PHONY : all vcfs tables reports
all : vcfs tables # reports

vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))

tables : $(TABLES)

reports : $(foreach type,gatk_indels gatk_snps,reports/$(type).dp_ft.grp)

filtered_snps : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.filtered.vcf)

vcf/%.gatk_snps.vcf : gatk/vcf/%.variants.snps.filtered.vcf
	$(INIT) ln -f $< $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.variants.indels.filtered.vcf
	$(INIT) ln -f $< $@

VARIANT_EVAL_GATK_REPORT = $(RSCRIPT) $(HOME)/share/scripts/variantEvalGatkReport.R

reports/%/index.html : reports/%.dp_ft.grp metrics/hs_metrics.txt
	$(call LSCRIPT,"$(VARIANT_EVAL_GATK_REPORT) --metrics $(word 2,$^) --outDir $(@D) $<")

include ~/share/modules/gatk.mk
