include modules/Makefile.inc

ifeq ($(VARIANT_MAPS_TEST),true)
include modules/test/variant_maps/project_config.inc
endif

LOGDIR = log/variant.maps.$(NOW)

# call R script
VARIANT_MAPS = $(RSCRIPT) modules/summary/variantMaps.R

# set input defaults
PROJECT_CONFIG ?= project_config.yaml
SAMPLES_CONFIG ?= samples.yaml
MUTS_TABLE_IN ?= summary/mutation_summary.xlsx
GENE_CN_IN ?= facets/geneCN.fill.txt
CNCF_PATH ?= facets/cncf
SEG_MAF_PATH ?= absolute/reviewed/SEG_MAF

# set run option defaults
VARIANT_MAP_INCLUDE_SILENT ?= true

variant_maps : summary/mutations.tsv summary/copy_number_alterations.tsv summary/variants.tsv

summary/mutations.tsv summary/copy_number_alterations.tsv summary/variants.tsv : \
project_config.yaml samples.yaml $(MUTS_TABLE_IN) $(GENE_CN_IN) $(CNCF_PATH) $(SEG_MAF_PATH) $(foreach pair,$(SAMPLE_PAIRS),$(CNCF_PATH)/$(pair).cncf.txt) $(foreach pair,$(SAMPLE_PAIRS),$(SEG_MAF_PATH)/$(pair)_ABS_MAF.txt)
	$(INIT) unset PYTHONPATH; source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); $(call LSCRIPT_CHECK_MEM,8G,30G, "$(VARIANT_MAPS) \
	--project_config=$(PROJECT_CONFIG) --samples_config=$(SAMPLES_CONFIG) --muts_table_in=$(MUTS_TABLE_IN) --gene_cn_in=$(GENE_CN_IN) --cncf_path=$(CNCF_PATH) --seg_maf_path=$(SEG_MAF_PATH) \
	--include_silent=$(VARIANT_MAP_INCLUDE_SILENT) \
	--muts_table_out=summary/mutations.tsv --cnas_table_out=summary/copy_number_alterations.tsv --variants_table_out=summary/variants.tsv")
