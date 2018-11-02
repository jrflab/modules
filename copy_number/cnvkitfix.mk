include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit

cnvkit : $(foreach sample,$(SAMPLES),cnvkit/$(sample).cnr)

define cnvkit-cnr
cnvkit/%.cnr : cnvkit/%.targetcoverage.cnn cnvkit/%.antitargetcoverage.cnn cnvkit/REFERENCE.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/REFERENCE.cnr -o cnvkit/$$(*).cnr")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))
				
.PHONY: $(PHONY)

