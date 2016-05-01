include modules/Makefile.inc

LOGDIR = log/variant.maps.$(NOW)

VARIANT_MAPS = $(RSCRIPT) modules/summary/variantMaps.R

# optional inclusion of silent mutations
ifneq ($(VARIANT_MAP_SILENT),)
VARIANT_MAPS_SILENT = --variant_maps_silent TRUE
endif

# define alternative mutation file
ifneq ($(VARIANT_MAP_MUT_FILE),)
VARIANT_MAPS_MUT_FILE = --variant_maps_silent TRUE
endif

summary/mutation_heatmap_variant.pdf summary/mutation_heatmap_variant_recurrent.pdf summary/mutation_heatmap_ccf.pdf summary/mutation_heatmap_ccf_recurrent.pdf : \
summary/mutation_summary.xlsx $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).cncf.txt) $(foreach pair,$(SAMPLE_PAIRS),absolute/reviewed/SEG_MAF/$(pair)_ABS_MAF.txt)
	$(INIT) unset PYTHONPATH; \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	$(call LSCRIPT_CHECK_MEM,8G,30G,"$(VARIANT_MAPS) --outFiles... $@ --inFiles... $<

