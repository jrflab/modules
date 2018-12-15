include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_fix_test.$(NOW)
PHONY += cnvkit cnvkit/cnr

cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).timestamp)

define cnvkit-cnr
cnvkit/cnr/%.timestamp : cnvkit/cnn/tumor/%.A.targetcoverage.cnn cnvkit/cnn/tumor/%.B.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"$(RSCRIPT) modules/test/copy_number/cnvkitfix.R --tumor $$(*) --normals '$$(NORMAL_SAMPLES)' && touch cnvkit/cnr/$$(*).timestamp")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))
				
.PHONY: $(PHONY)

