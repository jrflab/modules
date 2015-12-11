include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc


ALLTABLES_HIGH_MODERATE_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.low_modifier.novel.txt
ALLTABLES_SYNONYMOUS_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.synonymous.novel.txt
ALLTABLES_NONSYNONYMOUS_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.nonsynonymous.novel.txt
ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.low_modifier.novel.txt
ALLTABLES_SYNONYMOUS_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.synonymous.novel.txt
ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.nonsynonymous.novel.txt

ABSOLUTE_SOMATIC_TXTS ?= $(foreach sample,$(SAMPLE_SETS),absolute/tables/$(set).somatic.txt)
ABSOLUTE_SEGMENTS ?= $(foreach sample,$(SAMPLE_SETS),absolute/reviewed/SEG_MAF/$(set)_ABS_MAF.txt

mutation_summary: excel/mutation_summary.xlsx

EXCEL_MERGE_ABSOLUTE ?= false


ifeq ($(EXCEL_MERGE_ABSOLUTE),false)
excel/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_MUTECT) $(ALLTABLES_LOW_MODIFIER_MUTECT) $(ALLTABLES_SYNONYMOUS_MUTECT) $(ALLTABLES_NONSYNONYMOUS_MUTECT) $(ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN) $(ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN) $(ALLTABLES_SYNONYMOUS_STRELKA_VARSCAN) $(ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN)
	$(INIT) unset PYTHONPATH; \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/scripts/mutation_summary_excel.py $^ $@
else
excel/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_MUTECT) $(ALLTABLES_LOW_MODIFIER_MUTECT) $(ALLTABLES_SYNONYMOUS_MUTECT) $(ALLTABLES_NONSYNONYMOUS_MUTECT) $(ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN) $(ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN) $(ALLTABLES_SYNONYMOUS_STRELKA_VARSCAN) $(ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN) $(ABSOLUTE_SOMATIC_TXTS) $(ABSOLUTE_SEGMENTS)
	$(INIT) unset PYTHONPATH; \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/scripts/mutation_summary_excel.py --absolute_somatic_txts $(subst $(space),$(,),$(strip $(ABSOLUTE_SOMATIC_TXTS))) --absolute_segments $(subst $(space),$(,),$(strip $(ABSOLUTE_SEGMENTS))) $(wordlist 1,8,$^) $@
endif
