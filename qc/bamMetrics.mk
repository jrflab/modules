## defaults
LOGDIR = metrics/log

## includes
include modules/Makefile.inc
include modules/variant_callers/gatk.inc


.DELETE_ON_ERROR:

.SECONDARY: 

COLLECT_METRICS = $(JAVA) -Xmx8G -jar $(PICARD_DIR)/CollectMultipleMetrics.jar VALIDATION_STRINGENCY=LENIENT
COLLECT_WGS_METRICS = $(JAVA) -Xmx8G -jar $(PICARD_JAR) CollectWgsMetrics VALIDATION_STRINGENCY=LENIENT
COLLECT_GC_METRICS = $(JAVA) -Xmx8G -jar $(PICARD_DIR)/CollectGcBiasMetrics.jar VALIDATION_STRINGENCY=LENIENT

PHONY += bam_metrics
bam_metrics : summary_metrics gc flagstats wgs_metrics

PHONY += flagstats
flagstats : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats)
PHONY += summary_metrics
summary_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics)
PHONY += wgs_metrics
wgs_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics)
PHONY += dup
dup : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics)
PHONY += gc
gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics)

metrics/%.alignment_summary_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,13G,"$(COLLECT_METRICS) I=$< O=metrics/$* REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.wgs_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,13G,"$(COLLECT_WGS_METRICS) I=$< O=$@ REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.gc_bias_metrics : bam/%.bam
	$(call LSCRIPT_MEM,12G,13G,"$(COLLECT_GC_METRICS) I=$< O=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.flagstats : bam/%.bam
	$(call LSCRIPT_MEM,2G,3G,"$(SAMTOOLS) flagstat $< > $@")
	
bam/%.markdup.bam metrics/%.dup_metrics : bam/%.bam
	$(call LSCRIPT_MEM,18G,19G,"$(MARK_DUP) I=$< O=bam/$*.markdup.bam METRICS_FILE=metrics/$*.dup_metrics")

metrics/dup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics.txt)
	$(INIT) grep '^LIBRARY' $< > $@ && \
	for metrics in $^; do \
	    grep -A1 '^LIBRARY' $$metrics | sed '1d' >> $@; \
	done

.PHONY: $(PHONY)
