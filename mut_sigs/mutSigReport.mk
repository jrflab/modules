# annotate titan module
include modules/Makefile.inc
include modules/copy_number/titan.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/mutsig_report.$(NOW)

MUTSIG_REPORT = scripts/mutSigReport.Rmd
KNIT = $(RSCRIPT) scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/Alexandrov_NMF_signatures.txt

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : mutsig_report mutect_mutsig_report varscan_mutsig_report strelka_mutsig_report


mutsig_report : mutect_mutsig_report varscan_mutsig_report strelka_mutsig_report

mutect_mutsig_report : mutsig_report/mutect/index.html
varscan_mutsig_report : mutsig_report/varscan_snps/index.html
strelka_mutsig_report : mutsig_report/strelka_snps/index.html

define mutsig-report
mutsig_report/$1/index.html : $$(call VCFS,$1)
	$$(call LSCRIPT_MEM,8G,8G"$$(KNIT) $$(MUTSIG_REPORT) $$@ --name $1 --alexandrovData $$(ALEXANDROV_DATA) $$(if $$(TARGETS_FILE),--targetBed $$(TARGETS_FILE)) $$^")
endef
$(foreach type,varscan_snps strelka_snps mutect,$(eval $(call mutsig-report,$(type))))


