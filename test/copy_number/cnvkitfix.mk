include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit : $(foreach tumor_sample,$(TUMOR_SAMPLES),foreach normal_sample,$(NORMAL_SAMPLES),cnvkit/cnr/$(tumor_sample)_$(normal_sample).cnr)

define cnvkit-cnr
cnvkit/cnr/%.cnr : cnvkit/cnn/tumor/%.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn cnvkit/reference/%.cnr
	$$(call RUN,-c -s 6G -m 8G,"for i in $$(NORMAL_SAMPLES); do \
									cnvkit.py fix $$(<) $$(<<) 'cnvkit/reference/$$i.cnr' -o 'cnvkit/cnr/$$(*)_$$i.cnr'; \
								done && touch cnvkit/cnr/$$(*).timestamp")
	
endef
 $(foreach tumor_sample,$(TUMOR_SAMPLES),foreach normal_sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(tumor_sample),$(normal_sample))))
				
.PHONY: $(PHONY)
