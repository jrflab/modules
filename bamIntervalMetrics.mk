# generate bam interval metrics per sample

#NO_RM := true

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc
# picard format intervals file, needs requires sam format header

VPATH ?= bam

LOGDIR = log/metrics.$(NOW)

EXOME ?= false

ifeq ($(EXOME),true)
INTERVALS_FILE = $(HOME)/share/reference/SureSelect_50MB_S02972011_Regions_nochr_intervals.txt
else
INTERVALS_FILE ?= intervals.txt
endif

PLOT_HS_METRICS = $(RSCRIPT) $(HOME)/share/scripts/plotHsMetrics.R

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all hs_metrics amplicon_metrics report

all : hs_metrics report

hs_metrics : metrics/hs_metrics.txt metrics/interval_hs_metrics.txt

amplicon_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt)

report : metrics/interval_report/index.html

# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt : %.bam %.bam.bai
	$(call LSCRIPT_MEM,10G,20G,"$(CALC_HS_METRICS) INPUT=$< OUTPUT=$@ METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$(REF_FASTA) PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=$(INTERVALS_FILE) BAIT_SET_NAME=hs BAIT_INTERVALS=$(INTERVALS_FILE)")

# not sure how this differs from above, see picard doc
metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt : %.bam %.bam.bai
	$(call LSCRIPT_MEM,10G,20G,"$(COLLECT_TARGETED_METRICS) INPUT=$< REFERENCE_SEQUENCE=$(REF_FASTA) OUTPUT=$@ AMPLICON_INTERVALS=$(INTERVALS_FILE) TARGET_INTERVALS=$(INTERVALS_FILE) METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt")

# summarize metrics into one file
metrics/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1 > $@; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
	    sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics >> $@; \
	done

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


include ~/share/modules/processBamMD5.mk
