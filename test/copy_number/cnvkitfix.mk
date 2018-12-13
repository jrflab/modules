include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).timestamp)

define cnvkit-cnr
cnvkit/cnr/%.timestamp : cnvkit/cnn/tumor/%.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn $(wildcard cnvkit/reference/$(NORMAL_SAMPLES).cnr)
	$$(call RUN,-c -s 6G -m 8G,"for i in $$(NORMAL_SAMPLES); do \
									cnvkit.py fix $$(<) $$(<<) 'cnvkit/reference/$i.cnr' -o 'cnvkit/cnr/$$(*)_$i.cnr' && \
									touch cnvkit/cnr/$$(*).timestamp")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))
				
.PHONY: $(PHONY)
