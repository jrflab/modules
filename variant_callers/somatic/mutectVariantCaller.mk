# Run mutect on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

LOGDIR = log/mutect.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all mutect_vcfs mutect_tables ext_output mut_report

all : mutect_vcfs mutect_tables ext_output mut_report
mutect_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))
mutect_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt)) \
	$(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff).txt)
ext_output : $(foreach pair,$(SAMPLE_PAIRS),mutect/tables/$(pair).mutect.txt)
mut_report : mutect/report/report.timestamp mutect/lowAFreport/report.timestamp mutect/highAFreport/report.timestamp

include ~/share/modules/variant_callers/somatic/mutect.mk
