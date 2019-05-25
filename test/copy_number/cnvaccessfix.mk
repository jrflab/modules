include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_fix.$(NOW)
PHONY += cnvkit cnvaccess/cnr

cnvaccess_fix : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).A.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).B.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).C.cnr)

define cnvaccess-cnr
cnvaccess/cnr/%.A.cnr : cnvaccess/cnn/tumor/%.A.targetcoverage.cnn cnvaccess/cnn/tumor/%.A.antitargetcoverage.cnn cnvaccess/reference/on_target_pool_A.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/on_target_pool_A.cnr -o cnvaccess/cnr/$$(*).A.cnr")

cnvaccess/cnr/%.B.cnr : cnvaccess/cnn/tumor/%.B.targetcoverage.cnn cnvaccess/cnn/tumor/%.B.antitargetcoverage.cnn cnvaccess/reference/on_target_pool_B.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/on_target_pool_B.cnr -o cnvaccess/cnr/$$(*).B.cnr")
	
cnvaccess/cnr/%.C.cnr : cnvaccess/cnn/tumor/%.C.targetcoverage.cnn cnvaccess/cnn/tumor/%.C.antitargetcoverage.cnn cnvaccess/reference/off_target_no_pool.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/off_target_no_pool.cnr -o cnvaccess/cnr/$$(*).C.cnr")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-cnr,$(sample))))
				
.PHONY: $(PHONY)
