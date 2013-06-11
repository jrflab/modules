# Run somatic sniper on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/mutect.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

SOMATIC_SNIPER = /opt/common/somaticsniper/somaticsniper-1.0.2.2/bam-somaticsniper
SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

all : $(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).som_sniper.eff.nsfp.vcf)

define somsniper-tumor-normal
vcf/$1_$2.som_sniper.vcf : $1.bam $2.bam
	$$(call INIT_MEM,4G,8G) $$(SOMATIC_SNIPER) -F vcf -f $$(REF_FASTA) $$< $$(word 2,$$^) $$@
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call somsniper-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

# add pedigree info
# $(eval $(call pedigree-tumor-normal,tumor,normal))
define pedigree-tumor-normal
vcf/$1_$2.som_sniper.vcf : som_sniper/vcf/$1_$2.som_sniper.vcf
	$$(call INIT_MEM,6G,20G) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@; cat $< >> $$@
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call pedigree-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

include ~/share/modules/vcftools.mk
