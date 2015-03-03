# vim: set ft=make :
# Run Fastqc on bam files

include modules/Makefile.inc

FASTQC = $(PERL) /opt/common/FastQC/fastqc 
FASTQC_SUMMARY_PLOT = $(RSCRIPT) scripts/fastqcSummaryPlot.R

LOGDIR = log/fastqc.$(NOW)

.PHONY: all
.SECONDARY: 

all : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt) fastqc/all_summary.png

fastqc/%_fastqc/summary.txt : bam/%.bam
	$(call LSCRIPT_NAMED_MEM,$*_fastqc,4G,12G,"$(FASTQC) -o fastqc $^")

fastqc/all_summary.txt : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt)
	cut -f2 $< | tr '\n' '\t' | sed 's/^/Sample\t/; s/\t$$/\n/' > $@; \
	for sum in $^; do \
		sample=`head -1 $$sum | cut -f 3 | sed 's/\..*//'`; \
		cut -f1 $$sum | tr '\n' '\t' | sed "s/^/$$sample\t/; s/\t$$/\n/" >> $@; \
	done

fastqc/all_summary.png : fastqc/all_summary.txt
	$(FASTQC_SUMMARY_PLOT) --outFile $@ $<
