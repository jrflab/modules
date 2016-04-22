# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include modules/Makefile.inc

LOGDIR = log/exomeCNVLOH.$(NOW)
EXOMECNV = modules/copy_number/exomeCNV.R
EXOMECNVLOH = modules/copy_number/exomeCNVLOH.R
CREATE_BAF = $(PERL) $(HOME)/share/usr/bin/for.loh.files.pl

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all loh

all : loh

ifdef SAMPLE_PAIRS
LOH += $(foreach pair,$(SAMPLE_PAIRS),exomecnv/loh/$(pair).loh.txt)

define exomecnv-baf-tumor-normal-set
exomecnv/baf/$1_$2.baf_timestamp : vcf/$(subst $( ),_,$3).gatk_snps.target_ft.dp_ft.pass.vcf
	normal=`grep -m1 '^#CHROM' $$< | cut -f10- | tr '\t' '\n' | grep -n '^$2$$$$' | cut -f1 -d:`; \
	tumor=`grep -m1 '^#CHROM' $$< | cut -f10- | tr '\t' '\n' | grep -n '^$1$$$$' | cut -f1 -d:`; \
	$$(INIT) $$(CREATE_BAF) $$< exomecnv/baf/$1_$2.baf_1.txt exomecnv/baf/$1_$2.baf_2.txt $$$$normal $$$$tumor && touch $$@
exomecnv/baf/$1_$2.baf_1.txt : exomecnv/baf/$1_$2.baf_timestamp
exomecnv/baf/$1_$2.baf_2.txt : exomecnv/baf/$1_$2.baf_timestamp
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call exomecnv-baf-tumor-normal-set,$(tumor.$(pair)),$(normal.$(pair)),$(set.$(pair)))))


define exomecnv-loh-pair
exomecnv/loh/$1.loh.txt : exomecnv/baf/$1.baf_1.txt exomecnv/baf/$1.baf_2.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(RSCRIPT) $$(EXOMECNVLOH) --tumor $$< --normal $$(word 2,$$^) --outPrefix $$(@D)/$1")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call exomecnv-loh-pair,$(pair))))

else
LOH += $(foreach sample,$(SAMPLES),exomecnv/loh/$(sample).loh.txt)


endif

loh : $(LOH) 

include modules/vcf_tools/vcftools.mk
include modules/variant_callers/gatk.mk
