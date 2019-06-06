include modules/Makefile.inc

LOGDIR ?= log/fetchimpact.$(NOW)
PHONY += unprocessed_bam

fetch_impact : $(foreach sample,$(SAMPLES),unprocessed_bam/$(sample).bam)

define fetch-impact
unprocessed_bam/%.bam :
	$$(call RUN,-c -s 4G -m 12G,"scp luna.mskcc.org:/ifs/dmpshare/share/irb12_245/$$(*).bam unprocessed_bam/$$(*).bam")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call fetch-impact,$(sample))))

.PHONY : $(PHONY)
