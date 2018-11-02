include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_coverage.$(NOW)
PHONY += cnvkit

cnvkit : $(foreach sample,$(SAMPLES),cnvkit/$(sample).targetcoverage.cnn cnvkit/$(sample).antitargetcoverage.cnn)

define cnvkit-cnn
cnvkit/%.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) modules/bed/MSK-ACCESS-v1_0-probe-AB_cnvkit_ontarget.bed -o cnvkit/$$(*).targetcoverage.cnn")

cnvkit/%.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) modules/bed/MSK-ACCESS-v1_0-probe-AB_cnvkit_offtarget.bed -o cnvkit/$$(*).antitargetcoverage.cnn")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvkit-cnn,$(sample))))
		
.PHONY: $(PHONY)

