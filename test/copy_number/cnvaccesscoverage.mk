include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/bam cnvaccess/cnn cnvaccess/cnn/tumor cnvaccess/cnn/normal

cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/bam/$(sample).bam) \
					 #$(foreach sample,$(TUMOR_SAMPLES),cnvaccess/cnn/tumor/$(sample).pool-A.targetcoverage.cnn \
					 #cnvaccess/cnn/tumor/$(sample).pool-B.targetcoverage.cnn \
					 #cnvaccess/cnn/tumor/$(sample).no-pool.antitargetcoverage.cnn) \
					 #$(foreach sample,$(NORMAL_SAMPLES),cnvaccess/cnn/normal/$(sample).pool-A.targetcoverage.cnn \
					 #cnvaccess/cnn/normal/$(sample).pool-B.targetcoverage.cnn \
					 #cnvaccess/cnn/normal/$(sample).no-pool.antitargetcoverage.cnn)

ONTARGET_FILE_A = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
ONTARGET_FILE_B ?= ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFFTARGET_FILE = ~/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed

define mark-duplicates
cnvaccess/bam/%.bam : fgbio/%.merged.bam
	$$(call RUN,-c -n 1 -s 12G -m 18G,"java -Djava.io.tmpdir=$(TMPDIR) -Xmx16G -jar $$(PICARD_JAR) MarkDuplicates \
									   I=$$(<) \
									   O=cnvaccess/bam/$$(*).bam \
									   M=cnvaccess/bam/$$(*).txt && \
									   samtools index cnvaccess/bam/$$(*).bam && \
									   cp cnvaccess/bam/$$(*).bam.bai cnvaccess/bam/$$(*).bai")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call mark-duplicates,$(sample))))


define cnvaccess-tumor-cnn
cnvaccess/cnn/tumor/%.pool-A.targetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvaccess/cnn/tumor/$$(*).pool-A.targetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).pool-A.antitargetcoverage.cnn")
	
cnvaccess/cnn/tumor/%.pool-B.targetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvaccess/cnn/tumor/$$(*).pool-B.targetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).pool-B.antitargetcoverage.cnn")

cnvaccess/cnn/tumor/%.no-pool.antitargetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/tumor/$$(*).no-pool.antitargetcoverage.cnn && \
									 touch cnvaccess/cnn/tumor/$$(*).no-pool.targetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvaccess-tumor-cnn,$(sample))))
		
define cnvaccess-normal-cnn
cnvaccess/cnn/normal/%.pool-A.targetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvaccess/cnn/normal/$$(*).pool-A.targetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).pool-A.antitargetcoverage.cnn")
	
cnvaccess/cnn/normal/%.pool-B.targetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvaccess/cnn/normal/$$(*).pool-B.targetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).pool-B.antitargetcoverage.cnn")

cnvaccess/cnn/normal/%.no-pool.antitargetcoverage.cnn : cnvaccess/bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/normal/$$(*).no-pool.antitargetcoverage.cnn && \
									 touch cnvaccess/cnn/normal/$$(*).no-pool.targetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvaccess-normal-cnn,$(sample))))
		
.PHONY: $(PHONY)
