include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference.$(NOW)
PHONY += cnvkit

cnvkit : cnvkit/REFERENCE.cnr

cnvkit/REFERENCE.cnr : $(wildcard cnvkit/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 12G -m 16G,"cnvkit.py reference cnvkit/*coverage.cnn -f $(REF_FASTA) --no-edge -o cnvkit/REFERENCE.cnr")
		
.PHONY: $(PHONY)

