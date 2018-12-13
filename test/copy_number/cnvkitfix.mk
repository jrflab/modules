include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit : $(foreach pair,$(SAMPLE_PAIRS),cnvkit/cnr/$(pair).cnr)

define cnvkit-cnr
cnvkit/cnr/$1_$2.cnr : cnvkit/cnn/tumor/$1.targetcoverage.cnn cnvkit/cnn/tumor/$1.antitargetcoverage.cnn cnvkit/reference/$2.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/$2.cnr -o cnvkit/cnr/$1_$2.cnr")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call cnvkit-cnr,$(tumor.$(pair)),$(normal.$(pair)))))
				
.PHONY: $(PHONY)
