# Create tumour-normal dindel vcf files
include ~/share/modules/Makefile.inc

SAMPLE_PAIR_FILE ?= sample_pairs.txt

TUMOR_SAMPLES := $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES := $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
NSAMPLES = $(words $(TUMOR_SAMPLES))

VPATH = dindel/vcf

LOGDIR = log/dindel.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

VCFS = $(foreach sample,$(TUMOR_SAMPLES),dindel/vcf/$(sample).dindel.sorted.annotated.tnFiltered.vcf)

all : $(VCFS) $(addsuffix .idx,$(VCFS)))

include ~/share/modules/tnFilter.mk
include ~/share/modules/dindel.mk
