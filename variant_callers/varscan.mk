# Run VarScan
##### DEFAULTS ######


REF ?= hg19
LOGDIR = log/varscan.$(NOW)

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include modules/Makefile.inc

FIX_VARSCAN_VCF = $(PERL) scripts/fixVarscanVcf.pl
FP_FILTER = $(PERL) $(HOME)/share/usr/bin/fpfilter.pl
BAM_READCOUNT = $(HOME)/share/usr/bin/bam-readcount


VPATH ?= bam

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all vcfs copycalls segments cnv reports tables

VALIDATION ?= false

SNP_VCF_EFF_FIELDS += VAF
INDEL_VCF_EFF_FIELDS += VAF

ANN_TYPES = eff # annotated
EFF_TYPES = high_moderate low_modifier
VARIANT_TYPES = varscan_snps varscan_indels

MIN_MAP_QUAL ?= 1
VARSCAN_MIN_COVERAGE ?= 8
VARSCAN_MIN_READS2 ?= 2
VARSCAN_MIN_AVG_QUAL ?= 15
VARSCAN_MIN_VAR_FREQ ?= 0.03
VARSCAN_MIN_FREQ_FOR_HOM ?= 0.75
ifeq ($(VALIDATION),true)
VARSCAN_P_VALUE ?= 99e-02
VARSCAN_STRAND_FILTER ?= 1
else
VARSCAN_P_VALUE ?= 0.9
VARSCAN_STRAND_FILTER ?= 0
endif

override VARSCAN_OPTS = --min-coverage $(VARSCAN_MIN_COVERAGE) \
	--min-reads2 $(VARSCAN_MIN_READS2) \
	--min-avg-qual $(VARSCAN_MIN_AVG_QUAL) \
	--min-var-freq $(VARSCAN_MIN_VAR_FREQ) \
	--min-freq-for-hom $(VARSCAN_MIN_FREQ_FOR_HOM) \
	--p-value $(VARSCAN_P_VALUE) \
	--strand-filter $(VARSCAN_STRAND_FILTER) 

FILTER_SUFFIX := dp_ft.dgd_ft.encode_ft.nft
ifdef TARGETS_FILE
FILTER_SUFFIX := $(FILTER_SUFFIX).target_ft
endif
ANN_SUFFIX := pass.dbsnp.cosmic.nsfp.eff

ifeq ($(VALIDATION),true)
VCF_SUFFIX.varscan_snps := $(FILTER_SUFFIX).$(ANN_SUFFIX)
VCF_SUFFIX.varscan_indels := $(FILTER_SUFFIX).$(ANN_SUFFIX)
else
VCF_SUFFIX.varscan_snps := $(FILTER_SUFFIX).$(ANN_SUFFIX).chasm.fathmm
VCF_SUFFIX.varscan_indels := $(FILTER_SUFFIX).$(ANN_SUFFIX)
endif

ifeq ($(HRUN),true)
HRUN_FILTER ?= 1
VCF_SUFFIX.varscan_indels := $(VCF_SUFFIX.varscan_indels).hrun.hrun_ft
endif

VCF_SUFFIXES = $(foreach type,$(VARIANT_TYPES),$(type).$(VCF_SUFFIX.$(type)))
TABLE_SUFFIXES = $(foreach suff,$(VCF_SUFFIXES),$(suff).tab.novel $(suff).tab \
				 $(foreach eff,$(EFF_TYPES),\
				 $(suff).tab.$(eff).novel $(suff).tab.$(eff)))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(VCF_SUFFIXES),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(TABLE_SUFFIXES),tables/$(sample).$(suff).txt))
ALLTABLES = $(foreach suff,$(TABLE_SUFFIXES),alltables/all.$(suff).txt)

all : vcfs tables cnv
variants : vcfs tables
cnv : copycalls segments
vcfs : $(VCFS)
tables : $(TABLES) $(ALLTABLES)
reports : $(foreach type,$(VARIANT_TYPES),reports/$(type).$(FILTER_SUFFIX).grp)


ifeq ($(SPLIT_CHR),true)
define varscan-chr-type
varscan/chr_vcf/%.$1.$2.vcf : bam/%.bam bam/%.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) mpileup2$2 \
	<($$(SAMTOOLS) mpileup -r $1 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--output-vcf $$(VARSCAN_OPTS)  --vcf-sample-list $$* | $$(FIX_VARSCAN_VCF) -s $$* > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach type,snp indel,$(eval $(call varscan-chr-type,$(chr),$(type)))))

define merge-varscan-vcfs
varscan/vcf/$1.%.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1.$$(chr).%.vcf)
	$$(INIT) grep '^##' $$< > $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-varscan-vcfs,$(sample))))

else # no splitting by chr

define varscan-type
varscan/vcf/%.$1.vcf : bam/%.bam bam/%.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(VARSCAN) mpileup2$1 \
	<($$(SAMTOOLS) mpileup -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--output-vcf --vcf-sample-list $$* $$(VARSCAN_OPTS) | $$(FIX_VARSCAN_VCF) -s $$* > $$@")
endef
$(foreach type,snp indel, $(eval $(call varscan-type,$(chr),$(type))))
endif

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.vcf
	$(INIT) ln $< $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.vcf
	$(INIT) ln $< $@

include modules/variant_callers/gatk.mk
