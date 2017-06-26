# mutation signature report
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/mutsig_report.$(NOW)

VCF2VRANGES = $(RSCRIPT) modules/mut_sigs/vcf2VRanges.R
KNIT = $(RSCRIPT) modules/scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/sanger_30_mutsig_prob.txt
MUTSIG_REPORT = modules/mut_sigs/mutSigReport2.Rmd
MUTSIG_REPORT_OPTS = --name $(PROJECT_NAME) \
					 --alexandrovData $(ALEXANDROV_DATA) \
					 $(if $(TARGETS_FILE),--targetBed $(TARGETS_FILE))

SNV_TYPE ?= mutect

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mutect_mutsig_reports

mutect_mutsig_reports : mutsig_report/mutect/mutsig_report.timestamp

mutsig_report/mutect/mutsig_report.timestamp : $(foreach pair,$(SAMPLE_PAIRS),mutsig_report/vrange/$(pair).$(SNV_TYPE).ft.VRanges.Rdata)
	$(call LSCRIPT_ENV_NAMED_PARALLEL_MEM,$(MUTSIG_REPORT_ENV),mutect_mutsig_report,4,3G,5G,"$(KNIT) $(MUTSIG_REPORT) $(@D) --ncores 4 --outDir $(@D) $(MUTSIG_REPORT_OPTS) $^ && touch $@")

mutsig_report/vrange/%.VRanges.Rdata : vcf/%.vcf
	$(call LSCRIPT_ENV_MEM,$(MUTSIG_REPORT_ENV),7G,10G,"$(VCF2VRANGES) --genome $(REF) --outFile $@ $<")

include modules/vcf_tools/vcftools.mk
