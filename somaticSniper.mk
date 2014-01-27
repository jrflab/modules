# Run somatic sniper on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

SOMATIC_SNIPER = /opt/common/somaticsniper/somaticsniper-1.0.2.2/bam-somaticsniper
SOMATIC_SNIPER_OPTS ?= -q 1 -p
SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer

VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT DP DP4 VAQ BQ MQ AMQ SSC

LOGDIR = log/som_sniper.$(NOW)

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all somsniper_vcfs somsniper_tables

FILTER_SUFFIX := ss_dp_ft.pass.dbsnp.nsfp.chasm.fathmm.eff
EFF_TYPES = silent missense nonsilent_cds nonsilent
ANN_TYPES = eff # annotated
VCF_SUFFIXES = som_sniper.$(FILTER_SUFFIX)
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),som_sniper.$(FILTER_SUFFIX).tab.$(eff).pass.novel)

#VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).$(suff).vcf))
VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(suff).vcf))
all : somsniper_vcfs somsniper_tables
somsniper_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))
somsniper_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt)) \
	$(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff).txt)

define somsniper-tumor-normal
som_sniper/vcf/$1_$2.som_sniper.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,4G,8G,"$$(SOMATIC_SNIPER) $$(SOMATIC_SNIPER_OPTS) -F vcf -f $$(REF_FASTA) $$< $$(word 2,$$^) $$@")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call somsniper-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

# add pedigree info
# $(eval $(call pedigree-tumor-normal,tumor,normal))
define pedigree-tumor-normal
vcf/$1_$2.som_sniper.vcf : som_sniper/vcf/$1_$2.som_sniper.vcf
	$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call pedigree-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

include ~/share/modules/vcftools.mk
