## defaults
VPATH ?= bam
LOGDIR = log/rnaseq_metrics.$(NOW)

## includes
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

PLOT_RNASEQ_METRICS = $(RSCRIPT) modules/qc/plotRnaseqMetrics.R

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all report

COLLECT_RNASEQ_METRICS = $(JAVA) -Xmx7G -jar $(JARDIR)/CollectRnaSeqMetrics.jar VALIDATION_STRINGENCY=LENIENT
STRAND_SPECIFICITY ?= NONE

all : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics) metrics/all.rnaseq_metrics metrics/all.normalized_coverage.rnaseq_metrics report

report : metrics/rnaseq_report/index.html


metrics/%.rnaseq_metrics : bam/%.bam
	$(call RUN,-c -s 8G -m 12G,"$(COLLECT_RNASEQ_METRICS) REF_FLAT=$(GENE_REF_FLAT) RIBOSOMAL_INTERVALS=$(RIBOSOMAL_INTERVALS) STRAND_SPECIFICITY=$(STRAND_SPECIFICITY) INPUT=$< REFERENCE_SEQUENCE=$(REF_FASTA) OUTPUT=$@ CHART_OUTPUT=$@.pdf VERBOSITY=ERROR")

metrics/all.rnaseq_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics)
	grep '^PF_BASES' $< > $@ && for i in $^; do sample=`echo $$i | sed 's:.*/::; s/\..*//'`; grep -A1 '^PF_BASES' $$i | tail -1 | awk -v sample=$$sample 'BEGIN { OFS = "\t" } { $$23=sample; print $$0 }'  >> $@; done

metrics/all.normalized_coverage.rnaseq_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics)
	grep -A101 '^normalized_position' $< | cut -f1 > $@ && for i in $^; do sample=`echo $$i | sed 's:.*/::; s/\..*//'`; grep -A101 '^normalized_position' $$i | cut -f2 | sed "s/All_Reads/$$sample/" | paste $@ - > $@.tmp && mv $@.tmp $@; done

metrics/rnaseq_report/index.html : metrics/all.rnaseq_metrics metrics/all.normalized_coverage.rnaseq_metrics
	$(PLOT_RNASEQ_METRICS) --outDir $(@D) $^
