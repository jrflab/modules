include modules/Makefile.inc

LOGDIR = log/summary.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:

NO_HOTSPOTS ?= false

ifeq ($(NO_HOTSPOTS),false)
HOTSPOT_TABLE = alltables/allTN.hotspot.tab.txt
TABLES += $(HOTSPOT_TABLE)
endif

SUMMARY_CONFIG = $(wildcard summary_config.yaml)

SNV_TYPE ?= mutect
INDEL_TYPE ?= somatic_indels
SUMMARY_SNV_TYPE ?= $(word 1,$(SNV_TYPE))
SUMMARY_INDEL_TYPE ?= $(word 1,$(INDEL_TYPE))


SNV_TABLES = alltables/allTN.$(SUMMARY_SNV_TYPE).tab.high_moderate.txt alltables/allTN.$(SUMMARY_SNV_TYPE).tab.low_modifier.txt \
			 alltables/allTN.$(SUMMARY_SNV_TYPE).tab.synonymous.txt alltables/allTN.$(SUMMARY_SNV_TYPE).tab.nonsynonymous.txt
INDEL_TABLES = alltables/allTN.$(SUMMARY_INDEL_TYPE).tab.high_moderate.txt alltables/allTN.$(SUMMARY_INDEL_TYPE).tab.low_modifier.txt \
			   alltables/allTN.$(SUMMARY_INDEL_TYPE).tab.synonymous.txt alltables/allTN.$(SUMMARY_INDEL_TYPE).tab.nonsynonymous.txt
TABLES += $(SNV_TABLES) $(INDEL_TABLES)


# Add optional absolute results to excel
# the $(wildcard x) syntax is used to check for existence of file
MUTATION_SUMMARY_OPTS = $(if $(findstring false,$(EXOME)),--include_all) $(if $(HOTSPOT_TABLE), --hotspot $(HOTSPOT_TABLE)) $(if $(SUMMARY_CONFIG), --summary_config $(SUMMARY_CONFIG))
ABSOLUTE_SOMATIC_TXTS ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),absolute/tables/$(set).somatic.txt))
ABSOLUTE_SEGMENTS ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),absolute/reviewed/SEG_MAF/$(set)_ABS_MAF.txt))
ifneq ($(and $(ABSOLUTE_SOMATIC_TXTS),$(ABSOLUTE_SEGMENTS)),)
EXCEL_ABSOLUTE_PARAMS = --absolute_segments $(subst $(space),$(,),$(strip $(ABSOLUTE_SEGMENTS))) --absolute_somatic_txts $(subst $(space),$(,),$(strip $(ABSOLUTE_SOMATIC_TXTS)))
endif

# Add optional annotations file to excel, used to add any other desired columns
# to the excel file inner merge on TUMOR_SAMPLE CHROM POS REF ALT
EXCEL_ANNOTATION ?= $(wildcard summary/annotation.tsv)
ifneq ($(EXCEL_ANNOTATION),)
EXCEL_ANNOTATION_PARAMS = --annotation_tsv $(EXCEL_ANNOTATION)
endif

mutation_summary: summary/mutation_summary.xlsx

summary/mutation_summary.xlsx : $(SNV_TABLES) $(INDEL_TABLES) $(ABSOLUTE_SOMATIC_TXTS) $(ABSOLUTE_SEGMENTS) $(EXCEL_FACETS_LOH) $(EXCEL_ANNOTATION) $(if $(findstring false,$(NO_HOTSPOTS)),$(HOTSPOT_TABLE))
	$(INIT) python modules/summary/mutation_summary_excel.py $(MUTATION_SUMMARY_OPTS) --output_tsv_dir $(@D)/tsv $(EXCEL_ABSOLUTE_PARAMS) $(EXCEL_FACETS_LOH_PARAMS) $(EXCEL_ANNOTATION_PARAMS) $(SNV_TABLES) $(INDEL_TABLES) $@

include modules/vcf_tools/vcftools.mk
