include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference_test.$(NOW)
PHONY += cnvkit cnvkit/reference

cnvkit : $(foreach sample,$(NORMAL_SAMPLES),cnvkit/reference/$(sample).cnr)

define cnvkit-reference
cnvkit/reference/%.cnr : cnvkit/cnn/normal/%.targetcoverage.cnn cnvkit/cnn/normal/%.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py reference $$(<) $$(<<) --no-edge -o cnvkit/reference/$$(*).cnr")
	
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-reference,$(sample))))
		
.PHONY: $(PHONY)

