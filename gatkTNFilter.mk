# naive tumour-normal filter for gatk indels

include ~/share/modules/Makefile.inc

SAMPLE_PAIR_FILE = sample_pairs.txt

TUMOR_SAMPLES := $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES := $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
NSAMPLES = $(words $(TUMOR_SAMPLES))

VPATH = gatk/vcf

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

VCFS = $(foreach sample,$(TUMOR_SAMPLES),gatk/vcf/$(sample).indels.annotated.filtered.tnFiltered.vcf)
TABLES = $(foreach sample,$(TUMOR_SAMPLES),gatk/tables/$(sample).indels.annotated.filtered.tnFiltered.novel.txt)

all : $(VCFS) $(addsuffix .idx,$(VCFS)) $(TABLES)
	
include ~/share/modules/tnFilter.mk
include ~/share/modules/gatkVariantCaller.mk
