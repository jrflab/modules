include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/Makefile.inc


ALLTABLES_HIGH_MODERATE_MUTECT = alltables/allTN.$(call VCF_SUFFIXES,mutect).tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_MUTECT = alltables/allTN.$(call VCF_SUFFIXES,mutect).tab.low_modifier.novel.txt
ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.low_modifier.novel.txt

mutation_summary: excel/mutation_summary.xlsx


excel/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_MUTECT) $(ALLTABLES_LOW_MODIFIER_MUTECT) $(ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN) $(ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN)
	$(INIT) source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/scripts/mutation_summary_excel.py $< $(<<) $(<<<) $(<<<<) $@
