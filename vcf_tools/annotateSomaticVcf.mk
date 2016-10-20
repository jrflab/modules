# annotate mutect, strelka and varscan vcf files

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/annotate_somatic_vcf.$(NOW)

ifeq ($(CLUSTER_ENGINE),"PBS")
..DUMMY := $(shell cp $(SNVBOX_CONF) modules/external/SNVBox/snv_box.conf)
..DUMMY := $(shell python modules/scripts/launcher_sql_db.py modules/db/chasm_db.yaml)
..DUMMY := $(shell python modules/scripts/launcher_sql_db.py modules/db/fathmm_db.yaml)
..DUMMY := $(shell python modules/scripts/launcher_sql_db.py modules/db/ensembl-hs-core-85-37_db.yaml)
endif

SNV_TYPE ?= mutect_snps
INDEL_TYPE ?= somatic_indels
#strelka_indels varscan_indels strelka_varscan_indels
VARIANT_TYPES ?= $(SNV_TYPE) $(INDEL_TYPE)


DEPTH_FILTER ?= 5
HRUN ?= false
FFPE_NORMAL_FILTER ?= false

ANN_PATHOGEN ?= false
ANN_FACETS ?= false
ANN_MUT_TASTE ?= false
ifeq ($(ANN_PATHOGEN),true)
$(if $(or $(findstring b37,$(REF)),$(findstring hg19,$(REF))),,\
	$(error non-hg19/b37 pathogen annotation unsupported))
ANN_FACETS = true
ANN_MUT_TASTE = true
endif

SOMATIC_ANN1 = $(if $(findstring mm10,$(REF)),mgp_dbsnp,dbsnp) \
    eff $(if $(ANNOVAR_REF),$(ANNOVAR_REF)_multianno)\
    $(if $(findstring b37,$(REF)),cosmic gene_ann cn_reg clinvar exac_nontcga hotspot_ann)

ifeq ($(HRUN),true)
SOMATIC_INDEL_ANN2 += hrun
endif
ifeq ($(ANN_MUT_TASTE),true)
SOMATIC_INDEL_ANN2 += mut_taste
endif
SOMATIC_SNV_ANN2 = $(if $(findstring b37,$(REF)),nsfp chasm fathmm parssnp)

# indel/snv initial round of annotations
SOMATIC_ANN2 = $(if $(findstring indel,$1),$(SOMATIC_INDEL_ANN2),$(SOMATIC_SNV_ANN2))
ifeq ($(ANN_FACETS),true)
$(if $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)),,\
	$(error no facets data available for annotations))
SOMATIC_ANN2 += facets
endif

# apply depth filter to varscan and mutect
# fix vcf sample header for strelka
SOMATIC_INDEL_FILTER1 =
SOMATIC_SNV_FILTER1 = $(if $(findstring mutect,$1),som_ad_ft,\
    $(if $(findstring true,$(FFPE_NORMAL_FILTER)),ffpe_som_ad_ft))
# target filter
SOMATIC_FILTER1 = $(if $(TARGETS_FILE),target_ft) $(if $(findstring indel,$1),$(call SOMATIC_INDEL_FILTER1,$1),$(call SOMATIC_SNV_FILTER1,$1))

# filters run after initial round of annotations (but before final annotations)
SOMATIC_FILTER2 = cft common_ft
# hrun filter
SOMATIC_FILTER2 += $(if $(findstring indel,$1),\
            $(if $(findstring true,$(HRUN)),hrun_ft))

# final annotations (run last)
ifeq ($(ANN_PATHOGEN),true)
SOMATIC_INDEL_ANN3 += indel_pathogen
endif
SOMATIC_SNV_ANN3 = $(if $(findstring b37,$(REF)),snp_pathogen)

SOMATIC_ANN3 = $(if $(findstring indel,$1),$(SOMATIC_INDEL_ANN3),$(SOMATIC_SNV_ANN3))

PHONY += all somatic_vcfs merged_vcfs
all : somatic_vcfs merged_vcfs
merged_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf_ann/$(pair).merged.vcf.gz)
somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(type)_vcfs)

MERGE_VCF = $(PYTHON) modules/vcf_tools/merge_vcf.py
MERGE_SCRIPT = $(call LSCRIPT_MEM,6G,7G,"$(MERGE_VCF) --out_file $@ $^")
define somatic-merged-vcf
# first filter round
vcf/%.$1.ft.vcf : $$(if $$(strip $$(call SOMATIC_FILTER1,$1)),$$(foreach ft,$$(call SOMATIC_FILTER1,$1),vcf/%.$1.$$(ft).vcf),vcf/%.$1.vcf)
	$$(MERGE_SCRIPT)
# first annotation round
vcf/%.$1.ft.ann.vcf : $$(foreach ann,$$(call SOMATIC_ANN1,$1),vcf/%.$1.ft.$$(ann).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
vcf/%.$1.ft2.vcf : $$(foreach ft,$$(call SOMATIC_FILTER2,$1),vcf/%.$1.ft.ann.$$(ft).vcf)
	$$(MERGE_SCRIPT)
# post-filter after first annotation round
vcf/%.$1.ft2.ann2.vcf : $$(if $$(strip $$(call SOMATIC_ANN2,$1)),$$(foreach ann,$$(call SOMATIC_ANN2,$1),vcf/%.$1.ft2.$$(ann).vcf),vcf/%.$1.ft2.vcf)
	$$(MERGE_SCRIPT)
vcf_ann/%.$1.vcf : $$(if $$(strip $$(call SOMATIC_ANN3,$1)),$$(foreach ann,$$(call SOMATIC_ANN3,$1),vcf/%.$1.ft2.ann2.$$(ann).vcf),vcf/%.$1.ft2.ann2.vcf)
	$$(MERGE_SCRIPT)
PHONY += $1_vcfs
$1_vcfs : $$(foreach pair,$$(SAMPLE_PAIRS),vcf_ann/$$(pair).$1.vcf vcf/$$(pair).$1.ft2.ann2.vcf vcf/$$(pair).$1.ft2.vcf vcf/$$(pair).$1.ft.ann.vcf vcf/$$(pair).$1.ft.vcf)
endef
$(foreach type,$(VARIANT_TYPES),$(eval $(call somatic-merged-vcf,$(type))))

define somatic-merge-vcf-tumor-normal
vcf/$1_$2.%.reord.vcf.gz : vcf_ann/$1_$2.%.vcf.gz
	$$(call LSCRIPT_MEM,4G,5G,"$$(BCFTOOLS2) view -l 5 -s $1$$(,)$2 $$^ > $$@")
vcf_ann/$1_$2.merged.vcf.gz : $$(foreach type,$$(VARIANT_TYPES),vcf/$1_$2.$$(type).reord.vcf.gz vcf/$1_$2.$$(type).reord.vcf.gz.tbi)
	$$(call LSCRIPT_MEM,4G,6G,"$$(BCFTOOLS2) concat -O z -a -D $$(filter %.vcf.gz,$$^) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call somatic-merge-vcf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY) 

include modules/vcf_tools/vcftools.mk
