include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_coverage.$(NOW)
PHONY += cnvkit cnvkit/cnn cnvkit/cnn/tumor cnvkit/cnn/normal

cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn) $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn)

define cnvkit-tumor-cnn
cnvkit/cnn/tumor/%.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/tumor/$$(*).targetcoverage.cnn")

cnvkit/cnn/tumor/%.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/tumor/$$(*).antitargetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-tumor-cnn,$(sample))))
		
define cnvkit-normal-cnn
cnvkit/cnn/normal/%.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/normal/$$(*).targetcoverage.cnn")

cnvkit/cnn/normal/%.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/normal/$$(*).antitargetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-normal-cnn,$(sample))))
		
.PHONY: $(PHONY)

