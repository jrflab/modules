# vim: set ft=make :
# museq module for use by jsm.mk

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR = log/museq.$(NOW)

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
.PHONY : museq_vcfs museq_tables

all : museq_vcfs museq_tables

FILTER_SUFFIX := dp_ft.pass.dbsnp.nsfp.eff.chasm.fathmm.transfic
EFF_TYPES = silent missense nonsilent_cds nonsilent
VCF_SUFFIXES = museq.$(FILTER_SUFFIX).vcf
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),museq.$(FILTER_SUFFIX).tab.$(eff).novel.txt)

VCFS = $(foreach suff,$(VCF_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(suff)))
museq_vcfs : $(VCFS) $(addsuffix .idx,$(VCFS))
museq_tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff))) $(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.$(suff))

define museq-tumor-normal-chr
museq/chr_vcf/$1_$2.$3.museq.vcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call LSCRIPT_MEM,24G,60G,"$$(MUSEQ_CLASSIFY) tumour:$$< normal:$$(word 2,$$^) reference:$$(REF_FASTA) model:$$(MUSEQ_MODEL) --config $$(MUSEQ_CONFIG) --interval $3 --purity $$(TUMOR_PURITY) --out $$@ &> $$(LOG)")
endef
#$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call museq-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))
$(foreach chr,$(CHROMOSOMES), \
	$(foreach i,$(SETS_SEQ), \
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call museq-tumor-normal-chr,$(tumor),$(call get_normal,$(set.$i)),$(chr))))))


# merge museq chunks
define museq-tumor-normal
museq/vcf/$1_$2.museq.vcf : $$(foreach chr,$$(CHROMOSOMES),museq/chr_vcf/$1_$2.$$(chr).museq.vcf)
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ), \
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call museq-tumor-normal,$(tumor),$(call get_normal,$(set.$i)),$(chr)))))

vcf/%.museq.vcf : museq/vcf/%.museq.vcf
	$(INIT) $(FIX_MUSEQ_VCF) -R $(REF_FASTA) $< > $@ 2> $(LOG)

include ~/share/modules/vcf_tools/vcftools.mk
