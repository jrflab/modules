# Run TEQC R library on bams
# vim: set ft=make :
include modules/Makefile.inc

LOGDIR ?= teqc/log

.PHONY: teqc
.DELETE_ON_ERROR:
.SECONDARY:

teqc : teqc_report/index.html

teqc/%.Rdata : bam/%.bam bam/%.bam.bai
	$(call INIT_MEM,12G,14G) $(RSCRIPT) modules/qc/TEQC.R --ref=$(REF) --outFile $@ $< $(TARGETS_FILE) &> $(LOGDIR)/$(@F).log

teqc_report/index.html : $(foreach sample,$(SAMPLES),teqc/$(sample).Rdata)
	$(call INIT_MEM,12G,14G) $(MKDIR) teqc_report; $(RSCRIPT) modules/qc/TEQCreport.R --outDir=$(@D) $^

