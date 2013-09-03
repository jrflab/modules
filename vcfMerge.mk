# Merge vcf files

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all variant_eval gt_concordance


FILTER_SUFFIX := dp_ft
VARIANT_TYPES := mutect museq
MERGE_SUFFIX = $(subst $( ),_,$(VARIANT_TYPES))


all : $(foreach 


merged_vcf/%.$(MERGE_SUFFIX).vcf : vcf/%.$(FILTER_SUFFIX).vcf
	

