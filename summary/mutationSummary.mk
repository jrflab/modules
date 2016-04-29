include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/summary.$(NOW)


ALLTABLES_HIGH_MODERATE_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.high_moderate.txt
ALLTABLES_LOW_MODIFIER_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.low_modifier.txt
ALLTABLES_SYNONYMOUS_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.synonymous.txt
ALLTABLES_NONSYNONYMOUS_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.nonsynonymous.txt
ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.high_moderate.txt
ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.low_modifier.txt
ALLTABLES_SYNONYMOUS_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.synonymous.txt
ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.nonsynonymous.txt

# Add optional absolute results to excel
# the $(wildcard x) syntax is used to check for existence of file
ABSOLUTE_SOMATIC_TXTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/tables/$(set).somatic.txt))
ABSOLUTE_SEGMENTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/reviewed/SEG_MAF/$(set)_ABS_MAF.txt))
ifneq ($(ABSOLUTE_SOMATIC_TXTS)$(ABSOLUTE_SEGMENTS),)
EXCEL_ABSOLUTE_PARAMS = --absolute_segments $(subst $(space),$(,),$(strip $(ABSOLUTE_SEGMENTS))) --absolute_somatic_txts $(subst $(space),$(,),$(strip $(ABSOLUTE_SOMATIC_TXTS)))
endif

# Add optional annotations file to excel, used to add any other desired columns
# to the excel file inner merge on TUMOR_SAMPLE CHROM POS REF ALT
EXCEL_ANNOTATION ?= $(wildcard summary/annotation.tsv)
ifneq ($(EXCEL_ANNOTATION),)
EXCEL_ANNOTATION_PARAMS = --annotation_tsv $(EXCEL_ANNOTATION)
endif

# Set max for ExAC_AF column in MUTATION_SUMMARY sheet and tsv
EXCEL_MAX_EXAC_AF ?= 1

mutation_summary: summary/mutation_summary.xlsx

summary/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_MUTECT) $(ALLTABLES_LOW_MODIFIER_MUTECT) $(ALLTABLES_SYNONYMOUS_MUTECT) $(ALLTABLES_NONSYNONYMOUS_MUTECT) $(ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN) $(ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN) $(ALLTABLES_SYNONYMOUS_STRELKA_VARSCAN) $(ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN) $(ABSOLUTE_SOMATIC_TXTS) $(ABSOLUTE_SEGMENTS) $(EXCEL_FACETS_LOH) $(EXCEL_ANNOTATION)
	$(INIT) unset PYTHONPATH; \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/summary/mutation_summary_excel.py --max_exac_af $(EXCEL_MAX_EXAC_AF) --output_tsv_dir $(@D)/tsv $(EXCEL_ABSOLUTE_PARAMS) $(EXCEL_FACETS_LOH_PARAMS) $(EXCEL_ANNOTATION_PARAMS) $(wordlist 1,8,$^) $@
