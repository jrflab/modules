# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include ~/share/modules/Makefile.inc

REF ?= hg19
SAMPLE_PAIR_FILE ?= sample_pairs.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)

READ_LENGTH ?= 75

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

LOGDIR = log/exomeCNV.$(NOW)
EXOMECNV = $(HOME)/share/scripts/exomeCNV.R
EXOMECNVLOH = $(HOME)/share/scripts/exomeCNVLOH.R
CREATE_BAF = $(PERL) $(HOME)/share/usr/for.loh.files.pl

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all cnv loh

all : cnv loh
	
cnv : $(foreach tumor,$(TUMOR_SAMPLES),exomecnv/cnv/$(tumor)_$(normal_lookup.$(tumor)).cnv.txt)
loh : $(foreach tumor,$(TUMOR_SAMPLES),exomecnv/loh/$(tumor)_$(normal_lookup.$(tumor)).loh.txt)

metrics/%.read_len : %.bam
	$(INIT) $(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' > $@

define exomecnv-cnv-tumor-normal
exomecnv/cnv/$1_$2.cnv.txt : gatk/read_depth/$1.read_depth gatk/read_depth/$2.read_depth 
	$$(call INIT_PARALLEL_MEM,8,1G,1.5G) $$(RSCRIPT) $$(EXOMECNV) --numThreads 8 --readLen $$(READ_LENGTH) --outDir $$(@D) $$<.sample_interval_summary $$(word 2,$$^).sample_interval_summary &> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call exomecnv-cnv-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

define exomecnv-baf-tumor-normal
exomecnv/baf/$1.baf%txt exomecnv/baf/$2.baf%txt : vcf/$1_$2.mutect%vcf
	$(INIT) $(CREATE_BAF) $< exomecnv/baf/$2.baf.txt exomecnv/baf/$1.baf.txt 1 2
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call exomecnv-baf-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

define exomecnv-loh-tumor-normal
exomecnv/loh/$1_$2.loh.txt : exomecnv/baf/$1.baf.txt exomecnvbaf/$2.baf.txt
	$$(call INIT_MEM,4G,6G) $$(RSCRIPT) $$(EXOMECNVLOH) --outDir $$(@D) $$^ 
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call exomecnv-loh-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

include ~/share/modules/readDepth.mk
