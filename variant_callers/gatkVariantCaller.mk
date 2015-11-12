# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

##### DEFAULTS ######
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
EFF_TYPES = low_modifier high_moderate synonymous nonsynonymous
VARIANT_TYPES = gatk_snps gatk_indels

VALIDATION ?= false

FILTERS = dp_ft \
	$(if $(NORMAL_VCF),nft) \
	$(if $(TARGETS_FILE),target_ft) \
	pass eff \
	$(if $(findstring mm10,$(REF)),mgp_dbsnp,dbsnp) \
	$(if $(findstring b37,$(REF)),cosmic nsfp clinvar \
		$(if $(and $(findstring snps,$1),$(findstring false,$(VALIDATION))),chasm fathmm))

FILTER_SUFFIX = $1.$(subst $( ),.,$(strip $(FILTERS)))
VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(call FILTER_SUFFIX,$(type)))
TABLE_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(call FILTER_SUFFIX,$(type)).tab \
				 $(call FILTER_SUFFIX,$(type)).tab.novel \
				 $(foreach eff,$(EFF_TYPES),$(call FILTER_SUFFIX,$(type)).tab.$(eff).novel \
				 	$(call FILTER_SUFFIX,$(type)).tab.$(eff)))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(VCF_SUFFIXES),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(TABLE_SUFFIXES),tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(TABLE_SUFFIXES),alltables/all.$(suff).txt)


PHONY += gatk gatk_vcfs gatk_tables gatk_reports

gatk : gatk_vcfs gatk_tables # reports

gatk_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))

gatk_tables : $(TABLES)

gatk_reports : $(foreach type,gatk_indels gatk_snps,reports/$(type).dp_ft.grp)


include modules/variant_callers/gatk.mk

.PHONY : $(PHONY)
