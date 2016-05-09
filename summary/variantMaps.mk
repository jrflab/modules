include modules/Makefile.inc

LOGDIR = log/variant.maps.$(NOW)

# 
VARIANT_MAPS = $(RSCRIPT) modules/summary/variantMaps.R

# set input defaults
MUTS_TABLE_IN ?= summary/mutation_summary.xlsx
GENE_CN_IN ?= facets/geneCN.fill.txt
CNCF_PATH ?= facets/cncf
SEG_MAF_PATH ?= absolute/reviewed/SEG_MAF

# set run option defaults
VARIANT_MAP_INCLUDE_SILENT ?= true

mutations.tsv copy_number_alterations.tsv variants.tsv : \
project_config.yaml samples.yaml $(MUTS_TABLE_IN) facets/geneCN.fill.txt  $(foreach pair,$(SAMPLE_PAIRS),$(CNCF_PATH)/$(pair).cncf.txt) $(foreach pair,$(SAMPLE_PAIRS),$(SEG_MAF_PATH)/$(pair)_ABS_MAF.txt)
	$(INIT) unset PYTHONPATH; source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); $(call LSCRIPT_CHECK_MEM,8G,30G, "$(VARIANT_MAPS)
	--project_config=project_config.yaml --samples_config=samples.yaml --muts_table_in=$(MUTS_TABLE_IN) --gene_cn_in=$(GENE_CN_IN) --cncf_path=$(CNCF_PATH) --seg_maf_path=$(SEG_MAF_PATH) \
	--include_silent=$(VARIANT_MAP_INCLUDE_SILENT) \
	--muts_table_out=summary/mutations.tsv --cnas_table_out=summary/copy_number_alterations.tsv --variants_table_out=summary/variants.tsv")
