ifndef MUTECT_INC
VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT AD DP FA

#SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer -canon

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_OPTS = --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION)
MUTECT = $(JAVA) -Xmx11G -jar $(MUTECT_JAR) --analysis_type MuTect $(MUTECT_OPTS)

MUT_FREQ_REPORT = $(RSCRIPT) $(HOME)/share/scripts/plotSeqLogoFromMutect.R
endif
MUTECT_INC = true
