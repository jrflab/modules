# vim: set ft=make :
# Run Fastqc on bam files

include modules/Makefile.inc

FASTQC_SUMMARY_PLOT = $(RSCRIPT) modules/qc/fastqcSummaryPlot.R

LOGDIR ?= log/fastqc.$(NOW)

.PHONY: fastqc
.SECONDARY: 

fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt) fastqc/all_summary.txt

fastqc/%_fastqc.zip : bam/%.bam
	$(call LSCRIPT_NAMED_MEM,$*_fastqc,4G,12G,"$(FASTQC) -o fastqc $^")

fastqc/%_fastqc/summary.txt : fastqc/%_fastqc.zip
	$(INIT) $(UNZIP) -o -d fastqc $< &> $(LOG)

fastqc/all_summary.txt : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt)
	$(INIT) $(FASTQC_SUMMARY_PLOT) --outPrefix fastqc/all_summary $^ &> $(LOG)
