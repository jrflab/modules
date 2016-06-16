ifndef SOMATIC_VARIANT_CALLER_INC
DEPTH_FILTER ?= 5
HRUN ?= false
FFPE_NORMAL_FILTER ?= false
VALIDATION ?= false
ANN_FACETS ?= false

SOMATIC_ANN1 = $(if $(findstring mm10,$(REF)),mgp_dbsnp,dbsnp) \
    eff \
    $(if $(findstring b37,$(REF)),cosmic gene_ann cn_reg clinvar exac_nontcga hotspot_ann)

SOMATIC_INDEL_ANNS = $(if $(findstring b37,$(REF)),mut_taste) \
    $(if $(findstring true,$(HRUN)),hrun)
SOMATIC_SNV_ANNS = $(if $(findstring b37,$(REF)),nsfp chasm fathmm)
# pass filter for faster annotations
# indel/snv initial round of annotations
SOMATIC_ANN2 = $(if $(findstring indel,$1),$(SOMATIC_INDEL_ANNS),$(SOMATIC_SNV_ANNS))

# apply depth filter to varscan and mutect
# fix vcf sample header for strelka
SOMATIC_FILTER1 = $(if $(findstring varscan,$1)$(findstring mutect,$1),\
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft,som_ad_ft))
# target filter
SOMATIC_FILTER1 += $(if $(TARGETS_FILE),target_ft)

# filters run after initial round of annotations (but before final annotations)
SOMATIC_FILTER2 = cft common_ft
# hrun filter
SOMATIC_FILTER2 += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun_ft))

# final annotations (run last)
ifeq ($(ANN_FACETS),true)
SOMATIC_ANN2 += facets
ifeq ($(REF),b37)
SOMATIC_ANN3 += pathogen
endif
endif

MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_MEM,6G,7G,"$(MERGE_VCF) --out_file $@ $^")
define somatic-merged-vcf
# first filter round
vcf/ft/%.$1.ft.vcf : $$(foreach ft,$$(call SOMATIC_FILTER1,$1),vcf/ft/$$(ft)/%.$1.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# first annotation round
vcf/ann/%.$1.ann.vcf : $$(foreach ann,$$(call SOMATIC_ANN1,$1),vcf/ann/$$(ann)/%.$1.ft.$$(ann).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
vcf/ft2/%.$1.ft2.vcf : $$(foreach ft,$$(call SOMATIC_FILTER2,$1),vcf/ft2/$$(ft)/%.$1.ann.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
ifdef SOMATIC_ANN2
vcf/ann2/%.$1.ann2.vcf : $$(foreach ann,$$(call SOMATIC_ANN2,$1),vcf/ann2/$$(ann)/%.$1.ft2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf/ann2/%.$1.ann2.vcf : vcf/ft2/%.$1.ft2.pass.vcf
	$$(INIT) cp $$< $$@
endif
ifdef SOMATIC_ANN3
vcf_ann/%.$1.vcf : $$(foreach ann,$$(call SOMATIC_ANN3,$1),vcf/ann3/$$(ann)/%.$1.ann2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf_ann/%.$1.vcf: vcf/ann2/%.$1.ann2.vcf
	$$(INIT) cp $$< $$@
endif
endef
SOMATIC_VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).$1.vcf)
SOMATIC_MAFS = $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).$1.maf)

endif
SOMATIC_VARIANT_CALLER_INC = true
