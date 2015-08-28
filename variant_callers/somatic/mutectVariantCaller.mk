# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

LOGDIR = log/mutect.$(NOW)

PHONY += mutect mutect_vcfs mutect_tables ext_output mut_report

mutect : mutect_vcfs mutect_tables ext_output mut_report mutect_mutsig_report

include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/mutect.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/mut_sigs/mutSigReport.mk

..DUMMY := $(shell mkdir -p version; echo "$(MUTECT) &> version/mutect.txt")

mutect_vcfs : $(call VCFS,mutect) $(addsuffix .idx,$(call VCFS,mutect))
mutect_tables : $(call TABLES,mutect)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
mut_report : mutect/report/report.timestamp mutect/lowAFreport/report.timestamp mutect/highAFreport/report.timestamp
$(eval $(call mutsig-report-name-vcfs,mutect,$(call VCFS,mutect)))

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : $(PHONY)

