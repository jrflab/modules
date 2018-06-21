# generate bam interval metrics per sample

#NO_RM := true

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
# picard format intervals file, needs requires sam format header

VPATH ?= bam

LOGDIR ?= log/metrics.$(NOW)

PLOT_HS_METRICS = $(RSCRIPT) modules/qc/plotHsMetrics.R
NON_REF_FREQ = $(PERL) modules/qc/nonRefFreqFromPileup.pl
NON_REF_FREQ_BIN_SIZE = 0.01

SUMMARIZE_HS_METRICS = python modules/qc/summarize_hs_metrics.py
SUMMARIZE_IDXSTATS = python modules/qc/summarize_idxstats.py


.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: bam_interval_metrics hs_metrics amplicon_metrics interval_report non_ref_metrics insert_size_metrics idxstats

bam_interval_metrics : hs_metrics interval_report non_ref_metrics idxstats

non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.tsv)

hs_metrics : metrics/hs_metrics.tsv metrics/interval_hs_metrics.tsv metrics/hs_metrics.summary.tsv

amplicon_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.tsv)

interval_report : metrics/interval_report/interval_report.timestamp

insert_size_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).insert_size_metrics.tsv)

idxstats : metrics/idxstats_summary.tsv $(foreach sample,$(SAMPLES),metrics/$(sample).idxstats)

# interval metrics per sample
metrics/%.hs_metrics.tsv metrics/%.interval_hs_metrics.tsv : bam/%.bam bam/%.bam.bai
	$(call RUN,-s 10G -m 20G,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(CALC_HS_METRICS) INPUT=$< OUTPUT=metrics/$*.hs_metrics.tsv METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$(REF_FASTA) PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.tsv TARGET_INTERVALS=\$$TMP BAIT_SET_NAME=hs BAIT_INTERVALS=\$$TMP")

# not sure how this differs from above, see picard doc
metrics/%.amplicon_metrics.tsv metrics/%.interval_amplicon_metrics.tsv : bam/%.bam bam/%.bam.bai
	$(call RUN,-s 10G -m 20G,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\"  $(TARGETS_FILE) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(COLLECT_TARGETED_METRICS) INPUT=$< REFERENCE_SEQUENCE=$(REF_FASTA) OUTPUT=$@ AMPLICON_INTERVALS=\$$TMP TARGET_INTERVALS=\$$TMP METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.tsv")

# summarize interval metrics into one file
metrics/interval_hs_metrics.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.tsv)
	$(INIT) \
	sed '/^#/d; /^$$/d' $< | cut -f 1-6 > $@.tmp; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.interval_hs_metrics.tsv}); \
		sed '/^#/d; /^$$/d' $$metrics | cut -f 7,8 | sed "s/mean_coverage/$${samplename}_mean_coverage/; s/normalized_coverage/$${samplename}_normalized_coverage/" | paste $@.tmp - > $@; \
		cp $@ $@.tmp; \
	done; \
	rm -f $@.tmp

metrics/hs_metrics.summary.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.tsv)
	$(INIT) $(SUMMARIZE_HS_METRICS) --excel_file $(@:.tsv=.xlsx) --project_name $(PROJECT_NAME) $^ > $@ 2> $(LOG)

metrics/hs_metrics.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.tsv)
	$(INIT) \
		{ \
		sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1; \
		for metrics in $^; do \
			samplename=$$(basename $${metrics%%.hs_metrics.tsv}); \
			sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics; \
		done; \
		} > $@

metrics/interval_report/interval_report.timestamp : metrics/hs_metrics.tsv
	$(call RUN,-s 7G -m 10G,"$(PLOT_HS_METRICS) --outDir $(@D) $< && touch $@")

#metrics/%.interval_nonref_freq.tsv : bam/%.bam
#	$(call RUN,-s 8G -m 10G,"$(SAMTOOLS) mpileup -l $(TARGETS_FILE) -f $(REF_FASTA) $< | $(NON_REF_FREQ) -b $(NON_REF_FREQ_BIN_SIZE) > $@")

#metrics/%.insert_size_metrics.tsv : bam/%.bam
#	$(call RUN,-s 8G -m 10G,"$(call PICARD,CollectInsertSizeMetrics,8G) INPUT=$< OUTPUT=$@ \
#		REFERENCE_SEQUENCE=$(REF_FASTA) HISTOGRAM_FILE=$(@:.tsv=.pdf)")

#metrics/insert_size_metrics.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).insert_size_metrics.tsv)
#	$(INIT) \
#		{ \
#		sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/\s$$//; s/^/SAMPLE\t/' $< | head -1; \
#		for metrics in $^; do \
#			samplename=$$(basename $${metrics%%.insert_size_metrics.tsv}); \
#			grep -A1 '^MEDIAN_INSERT_SIZE' $$metrics | sed "1d; s/^/$$samplename\t/; s/\t\+$$//";  \
#		done; \
#		} > $@

#metrics/%.idxstats : bam/%.bam bam/%.bam.bai
#	$(call RUN,,"samtools idxstats $< > $@")

#metrics/idxstats_summary.tsv : $(foreach sample,$(SAMPLES),metrics/$(sample).idxstats)
#	$(INIT) $(SUMMARIZE_IDXSTATS) --excel_file $(@:.tsv=.xlsx) --project_name $(PROJECT_NAME) --targets_file $(TARGETS_FILE) $^ > $@ 2> $(LOG)

include modules/bam_tools/processBam.mk
