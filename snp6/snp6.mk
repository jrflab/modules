# SNP6 analysis module
# runs APT -> hapseg -> absolute

SNP6_USE_MAF ?= false

include modules/Makefile.inc

LOGDIR = log/snp6.$(NOW)

SNP6_LIB_DIR = $(HOME)/share/reference/SNP6_library_files
SNP6_CDF = $(SNP6_LIB_DIR)/GenomeWideSNP_6.Full.cdf
SNP6_MODELS = $(SNP6_LIB_DIR)/GenomeWideSNP_6.birdseed.models
SNP6_SPECIAL = $(SNP6_LIB_DIR)/GenomeWideSNP_6.Full.specialSNPs
SNP6_CHRX = $(SNP6_LIB_DIR)/GenomeWideSNP_6.chrXprobes
SNP6_CHRY = $(SNP6_LIB_DIR)/GenomeWideSNP_6.chrYprobes
SUMMARIZE_TARGET = $(HOME)/share/reference/penncnv_gw6/lib/hapmap.quant-norm.normalization-target.txt
APT_GENOTYPE = $(HOME)/share/usr/apt-1.15.0-x86_64-intel-linux/bin/apt-probeset-genotype
APT_SUMMARIZE = $(HOME)/share/usr/apt-1.15.0-x86_64-intel-linux/bin/apt-probeset-summarize

#SUMMARIZE_PATHWAY := $(shell echo "$(SUMMARIZE_PATHWAY)" | sed 's/\.[^,]\+//g; s/,/./g')
SUMMARIZE_PATHWAY := quant-norm.pm-only.med-polish.expr
# map summarize pathway to APT birdseed option
SUMMARIZE_PATHWAY_FULL.$(SUMMARIZE_PATHWAY) := quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true

GENOTYPE_PATHWAY = birdseed-v1

APT_SUMMARIZE_OPTS = --cdf-file $(SNP6_CDF) --target-sketch $(SUMMARIZE_TARGET)
APT_GENOTYPE_OPTS = -c $(SNP6_CDF) \
					--read-models-birdseed $(SNP6_MODELS) \
					--special-snps $(SNP6_SPECIAL) \
					--chrX-probes $(SNP6_CHRX) \
					--chrY-probes $(SNP6_CHRY) \
					--set-gender-method cn-probe-chrXY-ratio --write-models

HAPSEG = $(RSCRIPT) modules/snp6/hapseg.R
HAPSEG_PHASED_BGL_DIR = $(HOME)/share/reference/phasedBGL
DISEASE ?= breastcancer
HAPSEG_OPTS = --disease $(DISEASE) --phasedBGLDir $(HAPSEG_PHASED_BGL_DIR) --ref $(REF)
ABSOLUTE = $(RSCRIPT) modules/snp6/absolute.R
ABSOLUTE_OPTS = --disease $(DISEASE)

PENNCNV_AFFY = $(HOME)/share/usr/penncnv_gw6/bin/normalize_affy_geno_cluster.pl
PENNCNV_HAPMAP = $(HOME)/share/reference/penncnv_gw6/lib/hapmap.genocluster
PENNCNV_LOCFILE = $(HOME)/share/reference/penncnv_gw6/lib/affygw6.hg19.pfb

.SECONDARY:
.DELETE_ON_ERROR: 


ifdef SAMPLE_PAIRS
.PHONY: absoluteTN
absoluteTN : $(foreach pair,$(SAMPLE_PAIRS),absolute/$(pair).timestamp)

define hapseg-tumor-normal
hapseg/$1_$2/segdat.Rdata : apt/$$(GENOTYPE_PATHWAY).calls.txt apt/$$(GENOTYPE_PATHWAY).snp-models.txt apt/$$(SUMMARIZE_PATHWAY).summary.txt
	$$(call LSCRIPT_MEM,8G,10G,"$$(HAPSEG) $$(HAPSEG_OPTS) --callsFile $$(word 1,$$^) --clustersFile $$(word 2,$$^) --summaryFile $$(word 3,$$^) --resultsDir $$(@D) --outFile $$(@F) $1 $2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call hapseg-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define absolute-tumor-normal
absolute/$1_$2.timestamp : hapseg/$1_$2/segdat.Rdata
	$$(call LSCRIPT_MEM,8G,10G,"$$(ABSOLUTE) --tumour $1 --outPrefix segdat --resultsDir $$(@D) $$< && touch $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
else
.PHONY: absolute
absolute : $(foreach sample,$(SAMPLES),absolute/$(sample).timestamp)
endif

# APT birdseed-v1
apt/%.summary.txt : $(foreach sample,$(SAMPLES),cel/$(sample).CEL)
	$(call LSCRIPT_MEM,8G,10G,"$(APT_SUMMARIZE) -a $(SUMMARIZE_PATHWAY_FULL.$*) $(APT_SUMMARIZE_OPTS) -out-dir $(@D) $^")

apt/%.calls.txt apt/%.snp-models.txt :  $(foreach sample,$(SAMPLES),cel/$(sample).CEL)
	$(call LSCRIPT_MEM,8G,10G,"$(APT_GENOTYPE) $(APT_GENOTYPE_OPTS) -a $* --out-dir $(@D) $^")

hapseg/%/segdat.Rdata : apt/$(GENOTYPE_PATHWAY).calls.txt apt/$(GENOTYPE_PATHWAY).snp-models.txt apt/$(SUMMARIZE_PATHWAY).summary.txt
	$(call LSCRIPT_MEM,8G,10G,"$(HAPSEG) $(HAPSEG_OPTS) --callsFile $(word 1,$^) --clustersFile $(word 2,$^) --summaryFile $(word 3,$^) --resultsDir $(@D) --outFile $(@F) $*")

ifeq ($(SNP6_USE_MAF),true)
absolute/%.timestamp : hapseg/%/segdat.Rdata absolute/maf/%.maf.txt
	$(call LSCRIPT_MEM,8G,10G,"$(ABSOLUTE) --tumour $* --mafFile $(<<) --minMutAF 0 --outPrefix $* --resultsDir $(@D)/$* $< && touch $@")
else
absolute/%.timestamp : hapseg/%/segdat.Rdata
	$(call LSCRIPT_MEM,8G,10G,"$(ABSOLUTE) --tumour $* --outPrefix $* --resultsDir $(@D)/$* $< && touch $@")
endif
