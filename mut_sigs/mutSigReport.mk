# mutation signature report
include modules/Makefile.inc
include modules/variant_callers/variantCaller.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc


LOGDIR = log/mutsig_report.$(NOW)

MUTSIG_REPORT = modules/mut_sigs/mutSigReport.Rmd
KNIT = $(RSCRIPT) modules/scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/Alexandrov_NMF_signatures.txt

mutect_mutsig_reports : $(foreach pair,$(SAMPLE_PAIRS),mutsig_report/mutect/$(pair)_mutsig_report.timestamp)

mutsig_report/mutect/%_mutsig_report.timestamp : $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,mutect),vcf/%.$(suff).vcf)
	$(call LSCRIPT_NAMED_MEM,$*_mutect_mutsig_report,8G,40G,"$(KNIT) $(MUTSIG_REPORT) $(@D)/$* --outDir $(@D) --name $* --alexandrovData $(ALEXANDROV_DATA) $(if $(TARGETS_FILE),--targetBed $(TARGETS_FILE)) $^ && touch $@")

