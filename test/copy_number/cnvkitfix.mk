include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit : $(foreach sample,$(SAMPLES),cnvkit/cnr/$(sample).cnr)

define cnvkit-cnr
cnvkit/cnr/%.cnr : cnvkit/cnn/tumor/%.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn cnvkit/reference/combined_reference.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/combined_reference.cnr -o cnvkit/cnr/$$(*).cnr")
	
cnvkit/cnr/%.cnr : cnvkit/cnn/normal/%.targetcoverage.cnn cnvkit/cnn/normal/%.antitargetcoverage.cnn cnvkit/reference/combined_reference.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/combined_reference.cnr -o cnvkit/cnr/$$(*).cnr")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))
				
.PHONY: $(PHONY)

