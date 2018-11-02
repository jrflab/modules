include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference.$(NOW)
PHONY += cnvkit

TARGET_CNVKIT ?= $(wildcard cnvkit/$(NORMAL_SAMPLES).targetcoverage.cnn)
ANTITARGET_CNVKIT ?= $(wildcard cnvkit/$(NORMAL_SAMPLES).antitargetcoverage.cnn)

cnvkit : cnvkit/REFERENCE.cnr

cnvkit/REFERENCE.cnr : $(wildcard cnvkit/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 12G -m 16G,"cnvkit.py reference $(TARGET_CNVKIT) $(ANTITARGET_CNVKIT) -f $(REF_FASTA) --no-edge -o cnvkit/REFERENCE.cnr")
		
.PHONY: $(PHONY)

