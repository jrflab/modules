ifndef SOMATIC_VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN ?= false
FFPE_NORMAL_FILTER ?= false
VALIDATION ?= false

VALIDATION_FILTERS := dbsnp cosmic nsfp eff
INDEL_ANN_FILTERS := dbsnp cosmic nsfp eff
SNV_ANN_FILTERS := $(INDEL_ANN_FILTERS) chasm fathmm hap_insuf
# apply depth filter to varscan and mutect
FILTERS = $(if $(findstring varscan,$1) $(findstring mutect,$1),dp_ft \
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft))
# target filter
FILTERS += $(if $(TARGETS_FILE),target_ft)\
# hrun filter
FILTERS += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun hrun_ft))
# pass filter for faster annotations
FILTERS += pass
# indel/snv annotations
FILTERS += $(if $(findstring true,$(VALIDATION),$(VALIDATION_FILTERS),\
		   $(if $(findstring indel,$1),$(INDEL_ANN_FILTERS),$(SNV_ANN_FILTERS))))

FILTER_SUFFIX = $1.$(subst $( ),.,$(strip $(FILTERS)))

EFF_TYPES := silent missense nonsilent_cds nonsilent

VCF_SUFFIXES = $(foreach type,$1,$(call FILTER_SUFFIX,$(type)))

TABLE_EXTENSIONS = tab tab.novel $(foreach eff,$(EFF_TYPES),tab.$(eff) tab.$(eff).novel)

TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
				 $(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

VCFS = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call VCF_SUFFIXES,$1), \
        vcf/$(pair).$(suff).vcf))

TABLES = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call TABLE_SUFFIXES,$1),tables/$(pair).$(suff).txt))
TABLES += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/allTN.$(suff).txt)
endif
SOMATIC_VARIANT_CALLER_INC = true
