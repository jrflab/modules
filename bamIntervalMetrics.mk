# generate bam interval metrics per sample

#NO_RM := true

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc
# picard format intervals file, needs requires sam format header

VPATH ?= bam

LOGDIR = log/metrics.$(NOW)

EXOME ?= false

PLOT_HS_METRICS = $(RSCRIPT) $(HOME)/share/scripts/plotHsMetrics.R

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all hs_metrics amplicon_metrics report non_ref_metrics

all : hs_metrics report non_ref_metrics

non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.txt)

hs_metrics : metrics/hs_metrics.txt metrics/interval_hs_metrics.txt

amplicon_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt)

report : metrics/interval_report/index.html


# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt : %.bam %.bam.bai
	$(call LSCRIPT_MEM,10G,20G,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(CALC_HS_METRICS) INPUT=$< OUTPUT=$@ METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$(REF_FASTA) PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=$$TMP BAIT_SET_NAME=hs BAIT_INTERVALS=$$TMP")

# not sure how this differs from above, see picard doc
metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt : %.bam %.bam.bai
	$(call LSCRIPT_MEM,10G,20G,"TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\"  $(TARGETS_FILE) | awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(COLLECT_TARGETED_METRICS) INPUT=$< REFERENCE_SEQUENCE=$(REF_FASTA) OUTPUT=$@ AMPLICON_INTERVALS=$$TMP TARGET_INTERVALS=$$TMP METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt")

# summarize metrics into one file
metrics/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1 > $@; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
	    sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics >> $@; \
	done
# print mean, min, max for average bait coverage
# awk 'BEGIN {min = 10000000; max = 0; } NR > 1 {if ($21 > max) {max = $21; } if ($21 < min) {min = $21;} total += $21} END {print total / (NR-1), min, max }'

# summarize interval metrics into one file
metrics/interval_hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.txt)
	$(INIT) for metrics in $^; do \
		samplename=$$(basename $${metrics%%.interval_hs_metrics.txt}); \
		cut -f 7,8 $$metrics | sed "s/mean_coverage/$${samplename}_mean_coverage/; s/normalized_coverage/$${samplename}_normalized_coverage/" > $$metrics.tmp; \
	done; \
	cut -f 1-6 $< | paste - $(addsuffix .tmp,$^) > $@; \
	rm -f $(addsuffix .tmp,$^);


metrics/interval_report/index.html : metrics/hs_metrics.txt
	$(call LSCRIPT,"$(PLOT_HS_METRICS) --outDir $(@D) $<")

NON_REF_FREQ = $(PERL) $(HOME)/share/scripts/nonRefFreqFromPileup.pl
NON_REF_FREQ_BIN_SIZE = 0.01

metrics/%.interval_nonref_freq.txt : %.bam
	$(call LSCRIPT,"$(SAMTOOLS) mpileup -l $(TARGETS_FILE) -f $(REF_FASTA) $< | $(NON_REF_FREQ) -b $(NON_REF_FREQ_BIN_SIZE) > $@")


include ~/share/modules/processBamMD5.mk
