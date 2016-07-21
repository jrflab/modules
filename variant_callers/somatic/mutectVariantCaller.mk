# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

LOGDIR ?= log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_mafs ext_output mut_report

mutect : mutect_vcfs ext_output # mutect_mafs

include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/mutect.inc

..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

mutect_vcfs : $(foreach pair,$(SOMATIC_PAIRS),vcf/$(pair).mutect.vcf)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
mut_report : mutect/report/index.html mutect/lowAFreport/index.html mutect/highAFreport/index.html

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

