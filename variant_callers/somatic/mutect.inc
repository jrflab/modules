VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT AD DP FA

SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer -canon

MUTECT_JAR := $(HOME)/share/usr/lib/java/muTect-1.1.4.jar
MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
MUTECT = $(JAVA) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUT_FREQ_REPORT = $(RSCRIPT) $(HOME)/share/scripts/plotSeqLogoFromMutect.R


MUTECT_FILTER_SUFFIX := som_ad_ft
ifdef TARGETS_FILE
TARGET_FILTER ?= true
ifeq ($(TARGET_FILTER),true)
MUTECT_FILTER_SUFFIX := $(MUTECT_FILTER_SUFFIX).target_ft
endif
endif
MUTECT_FILTER_SUFFIX := $(MUTECT_FILTER_SUFFIX).pass.dbsnp.cosmic.nsfp.eff.chasm
EFF_TYPES = silent missense nonsilent_cds nonsilent
MUTECT_VCF_SUFFIXES = mutect.$(MUTECT_FILTER_SUFFIX)
MUTECT_TABLE_SUFFIXES = mutect.$(MUTECT_FILTER_SUFFIX).tab $(foreach eff,$(EFF_TYPES),mutect.$(MUTECT_FILTER_SUFFIX).tab.$(eff).novel mutect.$(MUTECT_FILTER_SUFFIX).tab.$(eff))
