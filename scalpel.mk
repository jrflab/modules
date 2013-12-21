# scalpel variant detection
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/scalpel.$(NOW)

SCALPEL = $(HOME)/share/usr/scalpel-0.1.1/scalpel
SCALPEL_OPTS = --ref $(REF_FASTA)
ifeq ($(EXOME),true)
SCALPEL_OPTS += --bed $(EXOME_BED)
endif
ifdef TARGETS_FILE
SCALPEL_OPTS += --bed $(TARGETS_FILE)
endif

SCALPEL2VCF = $(PERL) $(HOME)/share/scripts/scalpelToVcf.pl

FILTER_SUFFIX := dbsnp.nsfp.chasm.fathmm

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all


all : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).scalpel.$(FILTER_SUFFIX).vcf)

define scalpel-tumor-normal
scalpel/$1_$2/somatic.5x.indel.txt : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,8G,10G,"$$(SCALPEL) --somatic --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS) --dir $$(@D)")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define scalpel2vcf-tumor-normal
vcf/$1_$2.scalpel.vcf : scalpel/$1_$2/somatic.5x.indel.txt
	$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel2vcf-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

include ~/share/modules/vcftools.mk


