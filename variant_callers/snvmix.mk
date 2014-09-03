# call variants using SNVmix
include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = log/snvmix.$(NOW)

GET_PILEUP_RESULTS_BY_POS_SCRIPT = $(HOME)/share/scripts/getPileupResultsByPosns.pl
#CODON_FILE=/share/data/rmorin/data/codon_lookup.sort

THRESHOLD ?= 0.5
DEPTH_FILTER ?= 5
BASE_QUAL ?= 15
MAPPING_QUAL ?= 10

ARTIFACTFINDER_SCRIPT=${SCRIPTS_DIR}/generateArtifactFinder.pl
SNVMIX = $(HOME)/usr/bin/SNVMix2
SNVMIX_OPTS= -t MB -q ${BASE_QUAL} -Q ${MAPPING_QUAL}
SNVMIX_MODEL = $(HOME)/share/reference/shah_lobular_snvmix_model.txt

SNVMIX_TO_VCF = $(HOME)/share/scripts/snvmixToVCF.pl

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

FILTER_SUFFIX = dp_ft.dbsnp.nsfp
ifdef NORMAL_VCF
FILTER_SUFFIX := nft.$(FILTER_SUFFIX)
endif

ANN_TYPES = eff #annotated
EFF_TYPES = silent missense

VCF_SUFFIXES = $(foreach ann,$(ANN_TYPES),snvmix2.$(FILTER_SUFFIX).$(ann).vcf)
TABLE_SUFFIXES = $(foreach ann,$(ANN_TYPES),$(foreach eff,$(EFF_TYPES),snvmix2.$(FILTER_SUFFIX).$(ann).$(eff).pass.novel.txt))

VCFS = $(foreach sample,$(SAMPLES),$(foreach suff,$(VCF_SUFFIXES),vcf/$(sample).$(suff)))
TABLES = $(foreach sample,$(SAMPLES),$(foreach suff,$(TABLE_SUFFIXES),tables/$(sample).$(suff)))

all : $(VCFS) $(TABLES)

ifdef SNVMIX_MODEL
snvmix/tables/%.snvmix2.txt : bam/%.bam
	$(call LSCRIPT_MEM,8G,12G,"$(SAMTOOLS) mpileup -s -f $(REF_FASTA) $< | awk '\$$4 > 0' | $(SNVMIX) $(SNVMIX_OPTS) -m $(SNVMIX_MODEL) -p s -C  -o $@ 2> $(LOG)")
# TODO else train and classify
endif

vcf/%.snvmix2.vcf : snvmix/tables/%.snvmix2.txt
	$(INIT) $(PERL) $(SNVMIX_TO_VCF) -R $(REF_FASTA) -d $(DEPTH_FILTER) $< > $@ 2> $(LOG)

include ~/share/modules/vcf_tools/vcftools.mk
