include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit_coverage.$(NOW)
PHONY += cnvkit cnvkit/cnn cnvkit/cnn/tumor cnvkit/cnn/normal

cnvkit_coverage : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).A.targetcoverage.cnn cnvkit/cnn/tumor/$(sample).B.targetcoverage.cnn cnvkit/cnn/tumor/$(sample).C.targetcoverage.cnn) $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).A.targetcoverage.cnn cnvkit/cnn/normal/$(sample).B.targetcoverage.cnn cnvkit/cnn/normal/$(sample).C.targetcoverage.cnn)

define cnvkit-tumor-cnn
cnvkit/cnn/tumor/%.A.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvkit/cnn/tumor/$$(*).A.targetcoverage.cnn && \
									 touch cnvkit/cnn/tumor/$$(*).A.antitargetcoverage.cnn")
	
cnvkit/cnn/tumor/%.B.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvkit/cnn/tumor/$$(*).B.targetcoverage.cnn && \
									 touch cnvkit/cnn/tumor/$$(*).B.antitargetcoverage.cnn")

cnvkit/cnn/tumor/%.C.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/tumor/$$(*).C.antitargetcoverage.cnn && \
									 touch cnvkit/cnn/tumor/$$(*).C.targetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-tumor-cnn,$(sample))))
		
define cnvkit-normal-cnn
cnvkit/cnn/normal/%.A.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvkit/cnn/normal/$$(*).A.targetcoverage.cnn && \
									 touch cnvkit/cnn/normal/$$(*).A.antitargetcoverage.cnn")
	
cnvkit/cnn/normal/%.B.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvkit/cnn/normal/$$(*).B.targetcoverage.cnn && \
									 touch cnvkit/cnn/normal/$$(*).B.antitargetcoverage.cnn")

cnvkit/cnn/normal/%.C.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/normal/$$(*).C.antitargetcoverage.cnn && \
									 touch cnvkit/cnn/normal/$$(*).C.targetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-normal-cnn,$(sample))))
		
.PHONY: $(PHONY)

