ifndef SOMATIC_VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN ?= false
FFPE_NORMAL_FILTER ?= false
VALIDATION ?= false

SOMATIC_VALIDATION_ANNS := dbsnp eff \
    $(if $(findstring b37,$(REF)),cosmic nsfp)
SOMATIC_INDEL_ANNS := $(if $(findstring mm10,$(REF)),mgp_dbsnp,dbsnp) \
    eff \
    $(if $(findstring b37,$(REF)),cosmic nsfp gene_ann cn_reg clinvar exac_nontcga)
SOMATIC_SNV_ANNS := $(SOMATIC_INDEL_ANNS) \
    $(if $(findstring b37,$(REF)),chasm fathmm)
ifeq ($(ANN_FACETS),true)
SOMATIC_INDEL_ANNS += facets
SOMATIC_SNV_ANNS += facets
endif
# apply depth filter to varscan and mutect
# fix vcf sample header for strelka
SOMATIC_FILTERS = dp_ft \
    $(if $(findstring varscan,$1)$(findstring mutect,$1),\
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft))
# target filter
SOMATIC_FILTERS += $(if $(TARGETS_FILE),target_ft)
# hrun filter
SOMATIC_FILTERS += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun hrun_ft))
# pass filter for faster annotations
# indel/snv annotations
SOMATIC_ANNOTATIONS += $(if $(findstring true,$(VALIDATION)),$(SOMATIC_VALIDATION_ANNS),\
		   $(if $(findstring indel,$1),$(SOMATIC_INDEL_ANNS),$(SOMATIC_SNV_ANNS)))

SOMATIC_FILTER_SUFFIX = $1.$(subst $( ),.,$(strip $(call SOMATIC_FILTERS,$1)))
SOMATIC_ANN_SUFFIX = pass.$(subst $( ),.,$(strip $(call SOMATIC_ANNOTATIONS,$1)))
SOMATIC_SUFFIX = $(SOMATIC_FILTER_SUFFIX).$(SOMATIC_ANN_SUFFIX)

SOMATIC_EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous

SOMATIC_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATIC_SUFFIX,$(type)))

SOMATIC_TABLE_EXTENSIONS = tab tab.novel $(foreach eff,$(SOMATIC_EFF_TYPES),tab.$(eff) tab.$(eff).novel)

SOMATIC_TABLE_SUFFIXES = $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1),\
				 $(foreach ext,$(SOMATIC_TABLE_EXTENSIONS),$(suff).$(ext)))

SOMATIC_VCFS = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call SOMATIC_VCF_SUFFIXES,$1), \
        vcf/$(pair).$(suff).vcf))

SOMATIC_TABLES = $(foreach pair,$(SAMPLE_PAIRS),\
    $(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),tables/$(pair).$(suff).txt))
SOMATIC_TABLES += $(foreach suff,$(call SOMATIC_TABLE_SUFFIXES,$1),alltables/allTN.$(suff).txt)
endif
SOMATIC_VARIANT_CALLER_INC = true
