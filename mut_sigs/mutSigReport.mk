# mutation signature report
include modules/Makefile.inc
include modules/variant_callers/variantCaller.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc


LOGDIR = log/mutsig_report.$(NOW)

MUTSIG_REPORT = modules/mut_sigs/mutSigReport.Rmd
KNIT = $(RSCRIPT) modules/scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/sanger_30_mutsig_prob.txt

MAKE_TRINUC_MAF = $(PYTHON) $(HOME)/share/usr/mutation-signatures/make_trinuc_maf.py
MUT_SIG_MAIN = $(PYTHON) $(HOME)/share/usr/mutation-signatures/main.py
MUT_SIGS = $(HOME)/share/usr/mutation-signatures/Stratton_signatures30.txt

mutect_mutsig_reports : $(foreach pair,$(SAMPLE_PAIRS),mutsig_report2/mutect/$(pair)_mutsig_report.timestamp)

mutsig_report/mutect/%_mutsig_report.timestamp : $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,mutect),vcf/%.$(suff).vcf)
	$(call LSCRIPT_NAMED_MEM,$*_mutect_mutsig_report,8G,40G,"$(KNIT) $(MUTSIG_REPORT) $(@D)/$* --outDir $(@D) --name $* --alexandrovData $(ALEXANDROV_DATA) $(if $(TARGETS_FILE),--targetBed $(TARGETS_FILE)) $^ && touch $@")

mutsig_report2/mutect/%_mutsig_report.timestamp : $(foreach pair,maf/$(pair).$(call SOMATIC_FILTER_SUFFIX,mutect).pass.maf)
	$(call LSCRIPT_NAMED_MEM,$*_mutect_mutsig_report2,8G,12G,"$(MUT_SIG_MAIN) $(MUT_SIGS) $< $(@D)/$* && touch $@")
