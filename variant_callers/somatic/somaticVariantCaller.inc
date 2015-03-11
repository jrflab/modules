DEPTH_FILTER ?= 5

INDEL_ANN_FILTERS := dbsnp cosmic nsfp eff
SNV_ANN_FILTERS := $(INDEL_ANN_FILTERS) chasm fathmm hap_insuf
FILTERS = $(if $(DEPTH_FILTER),dp_ft \
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft)) \
    $(if $(TARGETS_FILE),target_ft)\
    $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN),hrun hrun_ft))) \
    pass
    $(if $(findstring indel,$1),$(INDEL_ANN_FILTERS)) \
    $(if $(findstring snv,$1),$(SNV_ANN_FILTERS))

FILTER_SUFFIX = $(subst $( ),.,$(strip $(FILTERS)))

EFF_TYPES := silent missense nonsilent_cds nonsilent

VCF_SUFFIXES = $(foreach type,$1,$(type).$(call FILTER_SUFFIX,$(type)))

TABLE_SUFFIXES = $(addsuffix .tab,$(VCF_SUFFIXES))
TABLE_SUFFIXES += $(foreach eff,$(EFF_TYPES),$(addsuffix .$(eff),$(TABLE_SUFFIXES)))
TABLE_SUFFIXES += $(addsuffix .novel,$(TABLE_SUFFIXES))
