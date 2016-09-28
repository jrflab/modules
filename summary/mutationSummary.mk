include modules/Makefile.inc

LOGDIR = log/summary.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:


SNV_TYPE ?= mutect_snps
INDEL_TYPE ?= somatic_indels

HOTSPOT_TABLE = alltables/allTN.hotspot.tab.txt
ALLTABLES_HIGH_MODERATE_SNVS = alltables/allTN.$(SNV_TYPE).tab.high_moderate.txt
ALLTABLES_LOW_MODIFIER_SNVS = alltables/allTN.$(SNV_TYPE).tab.low_modifier.txt
ALLTABLES_SYNONYMOUS_SNVS = alltables/allTN.$(SNV_TYPE).tab.synonymous.txt
ALLTABLES_NONSYNONYMOUS_SNVS = alltables/allTN.$(SNV_TYPE).tab.nonsynonymous.txt
ALLTABLES_HIGH_MODERATE_INDELS = alltables/allTN.$(INDEL_TYPE).tab.high_moderate.txt
ALLTABLES_LOW_MODIFIER_INDELS = alltables/allTN.$(INDEL_TYPE).tab.low_modifier.txt
ALLTABLES_SYNONYMOUS_INDELS = alltables/allTN.$(INDEL_TYPE).tab.synonymous.txt
ALLTABLES_NONSYNONYMOUS_INDELS = alltables/allTN.$(INDEL_TYPE).tab.nonsynonymous.txt

# Add optional absolute results to excel
# the $(wildcard x) syntax is used to check for existence of file
ABSOLUTE_SOMATIC_TXTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/tables/$(set).somatic.txt))
ABSOLUTE_SEGMENTS ?= $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/reviewed/SEG_MAF/$(set)_ABS_MAF.txt))
ifneq ($(and $(ABSOLUTE_SOMATIC_TXTS),$(ABSOLUTE_SEGMENTS)),)
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

summary/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_SNVS) $(ALLTABLES_LOW_MODIFIER_SNVS) $(ALLTABLES_SYNONYMOUS_SNVS) $(ALLTABLES_NONSYNONYMOUS_SNVS) $(ALLTABLES_HIGH_MODERATE_INDELS) $(ALLTABLES_LOW_MODIFIER_INDELS) $(ALLTABLES_SYNONYMOUS_INDELS) $(ALLTABLES_NONSYNONYMOUS_INDELS) $(HOTSPOT_TABLE) $(ABSOLUTE_SOMATIC_TXTS) $(ABSOLUTE_SEGMENTS) $(EXCEL_FACETS_LOH) $(EXCEL_ANNOTATION)
	$(INIT) python modules/summary/mutation_summary_excel.py --max_exac_af $(EXCEL_MAX_EXAC_AF) --output_tsv_dir $(@D)/tsv $(EXCEL_ABSOLUTE_PARAMS) $(EXCEL_FACETS_LOH_PARAMS) $(EXCEL_ANNOTATION_PARAMS) $(filter alltables/allTN.%,$^) $@

include modules/vcf_tools/vcftools.mk
