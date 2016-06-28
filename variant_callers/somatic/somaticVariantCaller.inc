ifndef SOMATIC_VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN_FILTER_THRESHOLD ?= 1
HRUN_FILTER ?= false
FFPE_NORMAL_FILTER ?= false
VALIDATION ?= false
ANN_FACETS ?= false

# apply depth filter to varscan and mutect
# fix vcf sample header for strelka
SOMATIC_FILTERS = $(if $(findstring tvc,$1),som_af_ft) \
$(if $(findstring varscan,$1)$(findstring mutect,$1),\
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft))
# target filter
SOMATIC_FILTERS += $(if $(TARGETS_FILE),target_ft)
SOMATIC_FILTER_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_FILTERS,$1)))

# first round of annotations after filtering
SOMATIC_ANNS = $(if $(findstring mm10,$(REF)),mgp_dbsnp,dbsnp) \
    eff $(if $(findstring b37,$(REF))$(findstring hg19,$(REF)),cosmic exac_nontcga) \
    $(if $(findstring indel,$1),\
        $(if $(findstring true,$(HRUN_FILTER)),hrun))
ifeq ($(ANN_FACETS),true)
SOMATIC_ANNS += facets
endif
SOMATIC_ANN_SUFFIX = $(subst $( ),.,$(strip pass $(call SOMATIC_ANNS,$1)))

# filters that are run after first round of annotations
SOMATIC_POST_FILTERS += cft common_ft
# hrun filter 
SOMATIC_POST_FILTERS += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN_FILTER)),hrun_ft))
SOMATIC_POST_FILTER_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_POST_FILTERS,$1) pass2))

SOMATIC_POST_ANNS += $(if $(findstring b37,$(REF))$(findstring hg19,$(REF)),\
    $(if $(findstring mutect,$1)$(findstring snp,$1),\
        nsfp chasm fathmm))
SOMATIC_POST_ANNS += $(if $(findstring b37,$(REF))$(findstring hg19,$(REF)),gene_ann cn_reg clinvar hotspot_ann)
ifeq ($(ANN_FACETS),true)
ifeq ($(REF),b37)
SOMATIC_POST_ANNS += $(if $(findstring strelka_indel,$1),mut_taste) pathogen
endif
endif
SOMATIC_POST_ANN_SUFFIX = $(subst $( ),.,$(strip $(call SOMATIC_POST_ANNS,$1)))

SOMATIC_SUFFIX = $(subst $( ),.,$(strip $1 $(call SOMATIC_FILTER_SUFFIX,$1) $(call SOMATIC_ANN_SUFFIX,$1) $(call SOMATIC_POST_FILTER_SUFFIX,$1) $(call SOMATIC_POST_ANN_SUFFIX,$1)))

SOMATIC_EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous

# vcfs w/o post filtering 
SOMATIC_FT1_SUFFIX = $(subst $( ),.,$(strip $1 $(SOMATIC_FILTER_SUFFIX) $(call SOMATIC_ANN_SUFFIX,$1)))
SOMATIC_FT1_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATIC_FT1_SUFFIX,$(type)))

SOMATIC_VCF_SUFFIXES = $(foreach type,$1,$(call SOMATIC_SUFFIX,$(type)))

SOMATIC_TABLE_EXTENSIONS = tab $(foreach eff,$(SOMATIC_EFF_TYPES),tab.$(eff))

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
