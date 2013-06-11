# GATK Somatic indel detector for tumour-normal matched samples
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

REF ?= hg19
LOGDIR = log/gatk_som_indel.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

VPATH = bam

LOGDIR = log/gatk.$(NOW)

INDEL_WINDOW_SIZE = 200

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).gatk_som_indels.eff.vcf)

som_indel_tables : $(foreach i,$(shell seq 1 $(NSAMPLES)),gatk/tables/$(word $i,$(NORMAL_SAMPLES))_$(word $i,$(TUMOR_SAMPLES)).som_indels.annotated.cbind.txt)

ifeq ($(SPLIT_CHR),true)
#$(call som-indel,tumor,normal,chr)
define som-indel-tumor-normal-chr
gatk/chr_vcf/$1_$2.$3.som_indels.vcf : gatk/chr_bam/$1.$3.realn.bam gatk/chr_bam/$2.$3.realn.bam gatk/chr_bam/$1.$3.realn.bai gatk/chr_bam/$2.$3.realn.bai
	$$(call INIT_MEM,8G,10G) $$(MKDIR) gatk/metrics; $$(call GATK_MEM,8G) -T SomaticIndelDetector -R $$(REF_FASTA) -I:tumor $$(word 1,$$^) -I:normal $$(word 2,$$^) -o $$@ --metrics_file gatk/metrics/$1_$2.som_indels.grp -L $3 --window_size $$(INDEL_WINDOW_SIZE) &> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call som-indel-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))

define merge-som-indel-tumor-normal
gatk/vcf/$1_$2.som_indels.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1_$2.$$(chr).som_indels.vcf)
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | vcfsorter.pl $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call merge-som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

else # no splitting by chromosome

#$(call som-indel-tumor-normal,tumor,normal)
define som-indel-tumor-normal
gatk/chr_vcf/$1_$2.som_indels.vcf : gatk/bam/$1.realn.bam gatk/bam/$2.realn.bam gatk/bam/$1.realn.bai gatk/bam/$2.realn.bai
	$$(call INIT_MEM,8G,10G) $$(MKDIR) gatk/metrics; $$(call GATK_MEM,8G) -T SomaticIndelDetector -R $$(REF_FASTA) -I:tumor $$(word 1,$$^) -I:normal $$(word 2,$$^) -o $$@ --metrics_file gatk/metrics/$1_$2.som_indels.grp --window_size $$(INDEL_WINDOW_SIZE) &> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
endif

#$(call pedigree-som-indel-tumor-normal,tumor,normal)
define pedigree-som-indel-tumor-normal
vcf/$1_$2.gatk_som_indels.vcf : gatk/vcf/$1_$2.som_indels.vcf
	$$(INIT) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@ && cat $$< >> $$@
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call pedigree-som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

include ~/share/modules/gatk.mk
