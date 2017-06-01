# annotate external vcfs
include modules/Makefile.inc

LOGDIR ?= log/ext_vcf.$(NOW)

EXT_NAME ?= ext

ext_ann : ext_vcfs ext_tables

ext_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).$(EXT_NAME).vcf)
ext_tables : alltables/allTN.$(EXT_NAME).tab.txt

LOGDIR ?= log/annotate_ext_vcf.$(NOW)

ANN_MUT_TASTE ?= false
ANN_PROVEAN ?= false
SOMATIC_ANN1 = fathmm chasm dbsnp hotspot_ann eff exac_nontcga cosmic clinvar cn_reg gene_ann nsfp $(ANNOVAR_REF)_multianno \
			   $(if $(findstring true,$(ANN_MUT_TASTE)),mut_taste) $(if $(findstring true,$(ANN_PROVEAN)),provean)
SOMATIC_ANN2 = $(if $(findstring true,$(ANN_PATHOGEN)),snp_pathogen indel_pathogen)

# target filter

PHONY += all vcfs
all : vcfs
vcfs : $(foreach type,$(VARIANT_TYPES),$(foreach sample,$(SAMPLES),vcf_ann/$(sample).$(type).vcf))

MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_CHECK_MEM,6G,7G,"$(MERGE_VCF) --ignore_filter $^ | $(VCF_SORT) $(REF_DICT) - > $@")

# first filter round
# first annotation round
vcf/%.$(EXT_NAME).ann.vcf : $(foreach ann,$(SOMATIC_ANN1),vcf/%.$(EXT_NAME).$(ann).vcf)
	$(MERGE_SCRIPT)
# second annotation round
vcf_ann/%.$(EXT_NAME).vcf : $(if $(strip $(SOMATIC_ANN2)),$(foreach ann,$(SOMATIC_ANN2),vcf/%.$(EXT_NAME).ann.$(ann).vcf),vcf/%.$(EXT_NAME).ann.vcf)
	$(MERGE_SCRIPT)


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

include modules/vcf_tools/vcftools.mk
include modules/variant_callers/gatk.inc
