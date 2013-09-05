# vim: set ft=make :
# museq module for use by jsm.mk

include ~/share/modules/Makefile.inc

LOGDIR = log/museq.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_FILE ?= samples.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

MUSEQ_PYTHON = $(HOME)/share/usr/anaconda/bin/python
MUSEQ_PYTHONPATH = $(HOME)/share/usr/anaconda/lib/python2.7/site-packages
MUSEQ_LD_LIBRARY_PATH = $(HOME)/share/usr/anaconda/lib

MUSEQ_TRAIN = PYTHONPATH=$(MUSEQ_PYTHONPATH) LD_LIBRARY_PATH=$(MUSEQ_LD_LIBRARY_PATH) $(MUSEQ_PYTHON) $(HOME)/share/usr/museq-3.0.0/train.py
MUSEQ_CLASSIFY = PYTHONPATH=$(MUSEQ_PYTHONPATH) LD_LIBRARY_PATH=$(MUSEQ_LD_LIBRARY_PATH) $(MUSEQ_PYTHON) $(HOME)/share/usr/museq-3.0.0/classify.py
MUSEQ_MODEL = $(HOME)/share/usr/museq-3.0.0/model.npz
MUSEQ_CONFIG = $(HOME)/share/usr/museq-3.0.0/metadata.config
TUMOR_PURITY = 70

FIX_MUSEQ_VCF = $(PERL) $(HOME)/share/scripts/fixMuseqVCF.pl

# VCF overrides for tumor-normal
VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT DP AD

SNP_EFF_FLAGS = -ud 0 -no-intron -no-intergenic -cancer -canon

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : museq_vcfs

all : museq_vcfs museq_tables

FILTER_SUFFIX := dp_ft.dbsnp.nsfp.fathmm.chasm
EFF_TYPES = silent missense nonsilent_cds nonsilent
ANN_TYPES = eff # annotated
VCF_SUFFIXES = $(foreach ann,$(ANN_TYPES),museq.$(FILTER_SUFFIX).$(ann).vcf)
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),$(foreach ann,$(ANN_TYPES),museq.$(FILTER_SUFFIX).$(ann).$(eff).pass.novel.txt))

museq_vcfs : $(foreach suff,$(VCF_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).$(suff)))
museq_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).$(suff))) $(foreach suff,$(TABLE_SUFFIXES),tables/allTN.$(suff))

define museq-tumor-normal-chr
museq/chr_vcf/$1_$2.$3.museq.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,8G,12G,"$$(MUSEQ_CLASSIFY) tumour:$$< normal:$$(word 2,$$^) reference:$$(REF_FASTA) model:$$(MUSEQ_MODEL) --config $$(MUSEQ_CONFIG) --interval $3 --purity $$(TUMOR_PURITY) --out $$@ &> $$(LOG)")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call museq-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))


# merge museq chunks
define museq-tumor-normal
museq/vcf/$1_$2.museq.vcf : $$(foreach chr,$$(CHROMOSOMES),museq/chr_vcf/$1_$2.$$(chr).museq.ann.vcf)
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call museq-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

vcf/%.museq.vcf : museq/vcf/%.museq.vcf
	$(INIT) $(FIX_MUSEQ_VCF) -R $(REF_FASTA) $< > $@ 2> $(LOG)

include ~/share/modules/vcftools.mk
