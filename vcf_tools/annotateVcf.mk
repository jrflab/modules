include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/annotate_vcf.$(NOW)
override VARIANT_TYPES = gatk_snps gatk_indels

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

PHONY += all vcfs
all : vcfs mafs
vcfs : $(foreach type,$(VARIANT_TYPES),$(type)_vcfs)


MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_MEM,6G,7G,"$(MERGE_VCF) --out_file $@ $^")
define merged-vcf
# first filter round
vcf/%.$1.ft.vcf : $$(foreach ft,$$(call FILTERS,$1),vcf/%.$1.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# first annotation round
vcf/%.$1.ann.vcf : $$(foreach ann,$$(call ANNOTATIONS,$1),vcf/%.$1.ft.$$(ann).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
ifdef POST_FILTERS
vcf/%.$1.ft2.vcf : $$(foreach ft,$$(call POST_FILTERS,$1),vcf/%.$1.ann.$$(ft).vcf)
	$$(MERGE_SCRIPT)
else
vcf/%.$1.ft2.vcf: vcf/%.$1.ann.vcf
	$$(INIT) cp $$< $$@
endif
# post-filter after first annotation round
ifdef POST_ANNS
vcf_ann/%.$1.vcf : $$(foreach ann,$$(call POST_ANNS,$1),vcf/%.$1.ft2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf_ann/%.$1.vcf: vcf/%.$1.ft2.pass.vcf
	$$(INIT) cp $$< $$@
endif
PHONY += $1_vcfs
$1_vcfs : $(foreach sample,$(SAMPLES),vcf_ann/$(sample).$1.vcf)
endef
$(foreach type,$(VARIANT_TYPES),$(eval $(call merged-vcf,$(type))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

include modules/vcf_tools/vcftools.mk
