# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include ~/share/modules/Makefile.inc

READ_LENGTH ?= 75


LOGDIR = log/exomeCNV.$(NOW)
OUTDIR ?= exomecnv
EXOMECNV = $(HOME)/share/scripts/exomeCNV.R
EXOMECNVLOH = $(HOME)/share/scripts/exomeCNVLOH.R
CREATE_BAF = $(PERL) $(HOME)/share/usr/bin/for.loh.files.pl

CBS_SENS_SPEC ?= 0.99
SENS_SPEC ?= 0.9999
ADMIXTURE_RATE ?= 0.5

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all cnv loh

all : cnv # loh
	
cnv : $(foreach pair,$(SAMPLE_PAIRS),exomecnv/cnv/$(pair).cnv.txt)
loh : $(foreach pair,$(SAMPLE_PAIRS),exomecnv/loh/$(pair).loh.txt)

metrics/%.read_len : %.bam
	$(INIT) $(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' > $@

define exomecnv-cnv-tumor-normal
$(OUTDIR)/cnv/$1_$2.cnv.txt : gatk/read_depth/$1.read_depth gatk/read_depth/$2.read_depth 
	$$(call LSCRIPT_MEM,8,1G,1.5G,"$$(RSCRIPT) $$(EXOMECNV) --cbsSensSpec $(CBS_SENS_SPEC) --sensSpec $(SENS_SPEC) --admixtureRate $(ADMIXTURE_RATE) --numThreads 8 --readLen $$(READ_LENGTH) --outDir $$(@D) $$<.sample_interval_summary $$(word 2,$$^).sample_interval_summary")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call exomecnv-cnv-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define exomecnv-baf-tumor-normal
$(OUTDIR)/baf/$1_$2.baf_timestamp : gatk/vcf/$1_$2.variants.snps.vcf
	$$(INIT) $$(CREATE_BAF) $$< $$(OUTDIR)/baf/$2.baf.txt $$(OUTDIR)/baf/$1.baf.txt 1 2
$$(OUTDIR)/baf/$1.baf.txt : $$(OUTDIR)/baf/$1_$2.baf_timestamp
$$(OUTDIR)/baf/$2.baf.txt : $$(OUTDIR)/baf/$1_$2.baf_timestamp
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call exomecnv-baf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


define exomecnv-loh-tumor-normal
$(OUTDIR)/loh/$1_$2.loh.txt : $(OUTDIR)/baf/$1.baf.txt $(OUTDIR)/baf/$2.baf.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(RSCRIPT) $$(EXOMECNVLOH) --outDir $$(@D) $$^") 
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call exomecnv-loh-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include ~/share/modules/qc/readDepth.mk
