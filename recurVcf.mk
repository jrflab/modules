# intersect vcf files
##### DEFAULTS ######

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

SET_VCF_SUFFIXES = gatk_snps.dp_ft.som_ft
PAIR_VCF_SUFFIXES = som_sniper.ss_dp_ft.ss_ft.pass mutect.som_ad_ft.pass
SAMPLE_SET_PAIR_VCF = $(foreach suff,$(SET_VCF_SUFFIXES),vcf/$(get_set.$1).$(suff).vcf) $(foreach suff,$(PAIR_VCF_SUFFIXES),vcf/$(get_pair.$1).$(suff).vcf)

RECUR_VCF = $(RSCRIPT) $(HOME)/share/scripts/recurVcf.R

all : $(foreach sample,$(TUMOR_SAMPLES),recur_pos/$(sample).recur.txt)

define recur-pos-tumor
recur_pos/$1.recur.txt : $$(call SAMPLE_SET_PAIR_VCF,$1)
	$$(INIT) $$(RECUR_VCF) --outFile $$@ $$^
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call recur-pos-tumor,$(tumor))))
