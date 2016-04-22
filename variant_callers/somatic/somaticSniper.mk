# Run somatic sniper on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

SOMATIC_SNIPER = /opt/common/somaticsniper/somaticsniper-1.0.2.2/bam-somaticsniper
SOMATIC_SNIPER_OPTS ?= -q 1 -p
SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer
FIX_AD = $(RSCRIPT) modules/variant_callers/somatic/somaticSniperFixAD.R

BAM_READCOUNT = $(HOME)/share/usr/bin/bam-readcount
FP_FILTER = $(PERL) $(HOME)/share/usr/bin/somsniper_fpfilter.pl

LOGDIR = log/som_sniper.$(NOW)

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all somsniper_vcfs somsniper_tables

ANN_SUFFIX = pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
FILTER_SUFFIX := ss_dp_ft.ss_ft.rn.som_ad_ft
ifdef TARGETS_FILE
FILTER_SUFFIX := $(FILTER_SUFFIX).target_ft
endif
EFF_TYPES = silent missense nonsilent_cds nonsilent
VCF_SUFFIX = som_sniper.$(FILTER_SUFFIX).$(ANN_SUFFIX)
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),$(VCF_SUFFIX).tab.$(eff).novel)

#VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).$(suff).vcf))
VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(VCF_SUFFIX).vcf)
all : somsniper_vcfs somsniper_tables
somsniper_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))
somsniper_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt)) \
	$(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff).txt)

define somsniper-tumor-normal
som_sniper/vcf/$1_$2.som_sniper.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,4G,8G,"$$(SOMATIC_SNIPER) $$(SOMATIC_SNIPER_OPTS) -F vcf -f $$(REF_FASTA) $$< $$(word 2,$$^) $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call somsniper-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

# add pedigree info
# $(eval $(call pedigree-tumor-normal,tumor,normal))
#$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
define pedigree-tumor-normal
vcf/$1_$2.som_sniper.vcf : som_sniper/vcf/$1_$2.som_sniper.fp.fixAD.vcf
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call pedigree-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

%.fixAD.vcf : %.vcf
	$(INIT) $(FIX_AD) --genome $(REF) --outFile $@ $<

define fpfilter-tumor-normal-chr
som_sniper/chr_vcf/$1_$2.$3.som_sniper.vcf : som_sniper/vcf/$1_$2.som_sniper.vcf
	$$(INIT) grep '^#' $$< > $$@ && grep -P '^$3\t' $$< >> $$@ || true

som_sniper/chr_vcf/$1_$2.$3.som_sniper.fp.vcf : som_sniper/chr_vcf/$1_$2.$3.som_sniper.vcf bam/$1.bam
	$$(call LSCRIPT_MEM,8G,35G,"$$(FP_FILTER) --output-basename $$@ --snp-file $$< --readcount-file <($$(BAM_READCOUNT) -f $$(REF_FASTA) $$(<<) $3) &> /dev/null && mv $$@.fp_pass $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach chr,$(CHROMOSOMES),\
		$(eval $(call fpfilter-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

som_sniper/vcf/%.som_sniper.fp.vcf : $(foreach chr,$(CHROMOSOMES),som_sniper/chr_vcf/%.$(chr).som_sniper.fp.vcf)
	$(INIT) grep '^#' $< > $@ && sed '/^#/d' $^ >> $@

include modules/vcf_tools/vcftools.mk
