include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference.$(NOW)
PHONY += cnvkit cnvkit/reference

cnvkit_reference : cnvkit/reference/combined_reference.cnr

cnvkit/reference/combined_reference.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvkit/cnn/normal/*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/combined_reference.cnr")
		
.PHONY: $(PHONY)

