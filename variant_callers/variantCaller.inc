ifndef VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN ?= false
VALIDATION ?= false

VALIDATION_ANNS := dbsnp eff \
    $(if $(findstring b37,$(REF)),cosmic nsfp)
INDEL_ANNS := $(if $(findstring GRCm38,$(REF)),mgp_dbsnp,dbsnp) \
    eff \
    $(if $(findstring b37,$(REF)),cosmic nsfp gene_ann cn_reg clinvar exac_nontcga)
SNV_ANNS := $(INDEL_ANN_FILTERS) \
    $(if $(findstring b37,$(REF)),chasm fathmm)

FILTERS = dp_ft
# filter out germline mutations if there's a NORMAL_VCF
FILTERS += $(if $(NORMAL_VCF),nft)
# target filter
FILTERS += $(if $(TARGETS_FILE),target_ft)
# hrun filter
FILTERS += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun hrun_ft))

# indel/snv annotations
ANNOTATIONS += $(if $(findstring true,$(VALIDATION)),$(VALIDATION_ANNS),\
		   $(if $(findstring indel,$1),$(INDEL_ANNS),$(SNV_ANNS)))

FILTER_SUFFIX = $1.$(subst $( ),.,$(strip $(call FILTERS,$1)))
ANN_SUFFIX = pass.$(subst $( ),.,$(strip $(call ANNOTATIONS,$1)))
SUFFIX = $(FILTER_SUFFIX).$(ANN_SUFFIX)

EFF_TYPES := high_moderate low_modifier synonymous nonsynonymous

VCF_SUFFIXES = $(foreach type,$1,$(call SUFFIX,$(type)))

TABLE_EXTENSIONS = tab tab.novel $(foreach eff,$(EFF_TYPES),tab.$(eff) tab.$(eff).novel)

TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
				 $(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(call VCF_SUFFIXES,$1),vcf/$(sample).$(suff).vcf))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(call TABLE_SUFFIXES,$1),tables/$(sample).$(suff).txt))
TABLES += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/all.$(suff).txt)

endif
VARIANT_CALLER_INC = true
