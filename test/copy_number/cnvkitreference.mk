include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_reference_test.$(NOW)
PHONY += cnvkit cnvkit/reference

cnvkit : $(foreach sample,$(NORMAL_SAMPLES),cnvkit/reference/$(sample).A.cnr cnvkit/reference/$(sample).B.cnr)

define cnvkit-reference
cnvkit/reference/%.A.cnr : cnvkit/cnn/normal/%.A.targetcoverage.cnn cnvkit/cnn/normal/%.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py reference cnvkit/cnn/normal/$$(*).A.targetcoverage.cnn cnvkit/cnn/normal/$$(*).antitargetcoverage.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/$$(*).A.cnr")
	
cnvkit/reference/%.B.cnr : cnvkit/cnn/normal/%.B.targetcoverage.cnn cnvkit/cnn/normal/%.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py reference cnvkit/cnn/normal/$$(*).B.targetcoverage.cnn cnvkit/cnn/normal/$$(*).antitargetcoverage.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/$$(*).B.cnr")
	
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-reference,$(sample))))
		
.PHONY: $(PHONY)

