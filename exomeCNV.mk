# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include ~/share/modules/Makefile.inc

REF ?= hg19
SAMPLE_PAIR_FILE ?= sample_pairs.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

LOGDIR = log/exomeCNV.$(NOW)
EXOMECNV = $(HOME)/share/scripts/exomeCNV.R

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach tumor,$(TUMOR_SAMPLES),exomecnv/$(tumor)_$(normal_lookup.$(tumor)).cnv.txt)

define exomecnv-tumor-normal
exomecnv/$1_$2.cnv.txt : gatk/read_depth/$1.read_depth gatk/read_depth/$2.read_depth
	$$(call INIT_PARALLEL_MEM,8,1G,1.5G) $$(RSCRIPT) $$(EXOMECNV) --numThreads 8 --outDir $$(@D) $$<.sample_interval_summary $$(word 2,$$^).sample_interval_summary &> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call exomecnv-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

include ~/share/modules/readDepth.mk
