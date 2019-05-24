include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit_fix : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).A.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).B.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).C.cnr)

define cnvkit-cnr
cnvkit/cnr/%.A.cnr : cnvkit/cnn/tumor/%.A.targetcoverage.cnn cnvkit/cnn/tumor/%.A.antitargetcoverage.cnn cnvkit/reference/on_target_pool_A.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/on_target_pool_A.cnr -o cnvkit/cnr/$$(*).A.cnr")

cnvkit/cnr/%.B.cnr : cnvkit/cnn/tumor/%.B.targetcoverage.cnn cnvkit/cnn/tumor/%.B.antitargetcoverage.cnn cnvkit/reference/on_target_pool_B.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/on_target_pool_B.cnr -o cnvkit/cnr/$$(*).B.cnr")
	
cnvkit/cnr/%.C.cnr : cnvkit/cnn/tumor/%.C.targetcoverage.cnn cnvkit/cnn/tumor/%.C.antitargetcoverage.cnn cnvkit/reference/off_target_no_pool.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/off_target_no_pool.cnr -o cnvkit/cnr/$$(*).C.cnr")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))
				
.PHONY: $(PHONY)

