include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cnn cnvaccess/cnn/tumor cnvaccess/cnn/normal

cnvaccess_coverage : $(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnn/tumor/$(sample).A.targetcoverage.cnn cnvaccess/cnn/tumor/$(sample).B.targetcoverage.cnn cnvaccess/cnn/tumor/$(sample).C.antitargetcoverage.cnn) $(foreach sample,$(NORMAL_SAMPLES),cnvaccess/cnn/normal/$(sample).A.targetcoverage.cnn cnvaccess/cnn/normal/$(sample).B.targetcoverage.cnn cnvaccess/cnn/normal/$(sample).C.antitargetcoverage.cnn)

define cnvaccess-tumor-cnn
cnvaccess/cnn/tumor/%.A.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvaccess/cnn/tumor/$$(*).A.targetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).A.antitargetcoverage.cnn")
	
cnvaccess/cnn/tumor/%.B.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvaccess/cnn/tumor/$$(*).B.targetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).B.antitargetcoverage.cnn")

cnvaccess/cnn/tumor/%.C.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/tumor/$$(*).C.antitargetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).C.targetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-tumor-cnn,$(sample))))
		
define cnvaccess-normal-cnn
cnvaccess/cnn/normal/%.A.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvaccess/cnn/normal/$$(*).A.targetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).A.antitargetcoverage.cnn")
	
cnvaccess/cnn/normal/%.B.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvaccess/cnn/normal/$$(*).B.targetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).B.antitargetcoverage.cnn")

cnvaccess/cnn/normal/%.C.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/normal/$$(*).C.antitargetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).C.targetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvaccess-normal-cnn,$(sample))))
		
.PHONY: $(PHONY)

