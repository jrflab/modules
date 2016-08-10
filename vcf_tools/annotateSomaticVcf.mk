# annotate mutect, strelka and varscan vcf files

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/annotate_somatic_vcf.$(NOW)

override VARIANT_TYPES = mutect strelka_indels varscan_indels strelka_varscan_indels

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

PHONY += all somatic_vcfs somatic_mafs
all : somatic_vcfs somatic_mafs
somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(type)_vcfs)

MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_MEM,6G,7G,"$(MERGE_VCF) --out_file $@ $^")
define somatic-merged-vcf
# first filter round
vcf/%.$1.ft.vcf : $$(foreach ft,$$(call SOMATIC_FILTER1,$1),vcf/%.$1.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# first annotation round
vcf/%.$1.ft.ann.vcf : $$(foreach ann,$$(call SOMATIC_ANN1,$1),vcf/%.$1.ft.$$(ann).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
vcf/%.$1.ft2.vcf : $$(foreach ft,$$(call SOMATIC_FILTER2,$1),vcf/%.$1.ft.ann.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
ifdef SOMATIC_ANN2
vcf/%.$1.ft2.ann2.vcf : $$(foreach ann,$$(call SOMATIC_ANN2,$1),vcf/%.$1.ft2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf/%.$1.ann2.vcf : vcf/%.$1.ft2.pass.vcf
	$$(INIT) cp $$< $$@
endif
ifdef SOMATIC_ANN3
vcf_ann/%.$1.vcf : $$(foreach ann,$$(call SOMATIC_ANN3,$1),vcf/%.$1.ft2.ann2.$$(ann).vcf)
	$$(MERGE_SCRIPT)
else
vcf_ann/%.$1.vcf: vcf/%.$1.ft2.ann2.vcf
	$$(INIT) cp $$< $$@
endif
PHONY += $1_vcfs
$1_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).$1.vcf)
endef
$(foreach type,$(VARIANT_TYPES),$(eval $(call somatic-merged-vcf,$(type))))

strelka_varscan_indels: $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).strelka_varscan_indels.vcf)
vcf_ann/%.strelka_varscan_indels.vcf : vcf_ann/%.varscan_indels.vcf vcf_ann/%.strelka_indels.vcf
	$(call LSCRIPT_MEM,9G,12G,"grep -P '^#' $< > $@ && $(BEDTOOLS) intersect -a $< -b $(<<) >> $@")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

include modules/vcf_tools/vcftools.mk
