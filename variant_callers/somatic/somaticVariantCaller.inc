DEPTH_FILTER ?= 5

INDEL_ANN_FILTERS := dbsnp cosmic nsfp eff
SNV_ANN_FILTERS := $(INDEL_ANN_FILTERS) chasm fathmm hap_insuf
FILTERS = $(if $(DEPTH_FILTER),dp_ft \
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft)) \
    $(if $(TARGETS_FILE),target_ft)\
    $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN),hrun hrun_ft))) \
    pass \
    $(if $(findstring indel,$1),$(INDEL_ANN_FILTERS)) \
    $(if $(findstring snp,$1),$(SNV_ANN_FILTERS))

FILTER_SUFFIX = $(subst $( ),.,$(strip $(FILTERS)))

EFF_TYPES := silent missense nonsilent_cds nonsilent

VCF_SUFFIXES = $(foreach type,$1,$(type).$(call FILTER_SUFFIX,$(type)))

TABLE_SUFFIXES = $(foreach suff,$(VCF_SUFFIXES),$(suff).tab $(suff).tab.novel \
    $(foreach eff,$(EFF_TYPES),$(suff).tab.$(eff) $(suff).tab.$(eff).novel))
