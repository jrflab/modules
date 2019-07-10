include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_fix.$(NOW)
PHONY += cnvkit cnvaccess/cnr

cnvaccess_fix : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).pool-A.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).pool-B.cnr) $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnr/$(sample).no-pool.cnr)

define cnvaccess-cnr
cnvaccess/cnr/%.pool-A.cnr : cnvaccess/cnn/tumor/%.pool-A.targetcoverage.cnn cnvaccess/cnn/tumor/%.pool-A.antitargetcoverage.cnn cnvaccess/reference/on_target_pool_A.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/on_target_pool_A.cnr -o cnvaccess/cnr/$$(*).pool-A.cnr")

cnvaccess/cnr/%.pool-B.cnr : cnvaccess/cnn/tumor/%.pool-B.targetcoverage.cnn cnvaccess/cnn/tumor/%.pool-B.antitargetcoverage.cnn cnvaccess/reference/on_target_pool_B.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/on_target_pool_B.cnr -o cnvaccess/cnr/$$(*).pool-B.cnr")
	
cnvaccess/cnr/%.no-pool.cnr : cnvaccess/cnn/tumor/%.no-pool.targetcoverage.cnn cnvaccess/cnn/tumor/%.no-pool.antitargetcoverage.cnn cnvaccess/reference/off_target_no_pool.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvaccess/reference/off_target_no_pool.cnr -o cnvaccess/cnr/$$(*).no-pool.cnr")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-cnr,$(sample))))
				
.PHONY: $(PHONY)
