include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/bam_metrics.$(NOW)
PHONY += metrics

COLLECT_METRICS = $(JAVA) -Xmx12G -jar $(PICARD_DIR)/CollectMultipleMetrics.jar VALIDATION_STRINGENCY=LENIENT
COLLECT_WGS_METRICS = $(JAVA) -Xmx12G -jar $(PICARD_JAR) CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT
COLLECT_GC_METRICS = $(JAVA) -Xmx12G -jar $(PICARD_DIR)/CollectGcBiasMetrics.jar VALIDATION_STRINGENCY=LENIENT

SUMMARIZE_IDXSTATS = python modules/qc/summarize_idxstats.py

bam_metrics : summary_metrics gc flagstats wgs_metrics

PHONY += flagstats
flagstats : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats)
PHONY += summary_metrics
summary_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics)
PHONY += wgs_metrics
wgs_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics) metrics/wgs_metrics_summary.tsv
PHONY += dup
dup : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics)
PHONY += gc
gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics)

metrics/%.alignment_summary_metrics : bam/%.bam
	$(call RUN,-s 18G -m 24G -w 7200,"$(COLLECT_METRICS) I=$< O=metrics/$(*).alignment_summary_metrics REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/wgs_metrics_summary.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics)
	$(INIT) (grep GENOME_TERRITORY $< | sed 's/^/SAMPLE\t/'; for x in $(SAMPLES); do grep -A1 GENOME_TERRITORY metrics/$$x.wgs_metrics | sed 1d | sed "s/^/$$x\t/" ; done) > $@

metrics/%.wgs_metrics : bam/%.bam
	$(call RUN,-s 18G -m 24G -w 7200,"$(COLLECT_WGS_METRICS) I=$< O=$@ REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.gc_bias_metrics : bam/%.bam
	$(call RUN,-s 18G -m 24G -w 7200,"$(COLLECT_GC_METRICS) I=$< O=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.flagstats : bam/%.bam
	$(call RUN,-s 18G -m 24G -w 7200,"$(SAMTOOLS) flagstat $< > $@")
	
bam/%.markdup.bam metrics/%.dup_metrics : bam/%.bam
	$(call RUN,-s 18G -m 24G -w 7200,"$(MARK_DUP) I=$< O=bam/$*.markdup.bam METRICS_FILE=metrics/$*.dup_metrics")

metrics/dup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics.txt)
	$(INIT) grep '^LIBRARY' $< > $@ && \
	for metrics in $^; do \
	    grep -A1 '^LIBRARY' $$metrics | sed '1d' >> $@; \
	done

.PHONY: $(PHONY)
