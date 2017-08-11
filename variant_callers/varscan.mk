# Run VarScan
##### DEFAULTS ######

LOGDIR = log/varscan.$(NOW)

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include modules/Makefile.inc

FIX_VARSCAN_VCF = $(PERL) modules/variant_callers/fixVarscanVcf.pl
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
EFF_TYPES = high_moderate low_modifier synonymous nonsynonymous
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

VARIANT_TYPES = varscan_snps varscan_indels
VCFS = $(foreach sample,$(SAMPLES),$(foreach type,$(VARIANT_TYPES),vcf/$(sample).$(type).vcf))

all : vcfs
vcfs : $(VCFS)


ifeq ($(SPLIT_CHR),true)
define varscan-chr-type
varscan/chr_vcf/%.$1.$2.vcf : bam/%.bam bam/%.bam.bai
	$$(call RUN,-s 9G -m 12G,"$$(VARSCAN) mpileup2$2 \
	<($$(SAMTOOLS) mpileup -r $1 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
	--output-vcf $$(VARSCAN_OPTS)  --vcf-sample-list $$* | $$(FIX_VARSCAN_VCF) -s $$* > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach type,snp indel,$(eval $(call varscan-chr-type,$(chr),$(type)))))

define merge-varscan-vcfs
varscan/vcf/$1.%.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1.$$(chr).%.vcf)
	$$(call RUN,-s 4G -m 5G,"grep '^##' $$< > $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-varscan-vcfs,$(sample))))

else # no splitting by chr

define varscan-type
varscan/vcf/%.$1.vcf : bam/%.bam bam/%.bam.bai
	$$(call RUN,-s 9G -m 12G,"$$(VARSCAN) mpileup2$1 \
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
