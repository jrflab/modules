ifndef VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN ?= false
VALIDATION ?= false

ANNOTATIONS := $(if $(findstring GRCm38,$(REF)),mgp_dbsnp,dbsnp) \
    eff \
    $(if $(findstring b37,$(REF)),cosmic gene_ann cn_reg clinvar exac_nontcga)
#SNV_ANNS := $(INDEL_ANNS) \
#    $(if $(findstring b37,$(REF)),chasm fathmm)

FILTERS = dp_ft
# filter out germline mutations if there's a NORMAL_VCF
FILTERS += $(if $(NORMAL_VCF),nft)
# target filter
FILTERS += $(if $(TARGETS_FILE),target_ft)
# hrun filter
FILTERS += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun hrun_ft))

#POST_FILTERS += cft common_ft



MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_MEM,6G,7G,"$(MERGE_VCF) --out_file $@ $^")
define merged-vcf
# first filter round
vcf/ft/%.$1.ft.vcf : $$(foreach ft,$$(call FILTERS,$1),vcf/ft/$$(ft)/%.$1.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# first annotation round
vcf/ann/%.$1.ann.vcf : $$(foreach ann,$$(call ANNOTATIONS,$1),vcf/ann/$$(ann)/%.$1.ft.$$(ann).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
ifdef POST_FILTERS
vcf/ft2/%.$1.ft2.vcf : $$(foreach ft,$$(call POST_FILTERS,$1),vcf/ft2/$$(ft)/%.$1.ann.$$(ft).vcf)
	$$(MERGE_SCRIPT)
else
vcf/ft2/%.$1.ft2.vcf: vcf/ann/%.$1.ann.vcf
	$$(INIT) cp $$< $$@
endif
# post-filter after first annotation round
ifdef POST_ANNS
vcf_ann/%.$1.vcf : $$(foreach ann,$$(call POST_ANNS,$1),vcf/ann2/$$(ann)/%.$1.ft2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf_ann/%.$1.vcf: vcf/ft2/%.$1.ft2.pass.vcf
	$$(INIT) cp $$< $$@
endif
endef
VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).$1.vcf)
MAFS = $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).$1.maf)

endif
VARIANT_CALLER_INC = true
