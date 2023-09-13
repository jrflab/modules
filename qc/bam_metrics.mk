include modules/Makefile.inc

LOGDIR ?= log/bam_metrics.$(NOW)

bam_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).idx_stats.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).aln_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).insert_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).gc_metrics.txt) \
	      summary/idx_metrics.txt \
	      summary/aln_metrics.txt \
	      summary/insert_metrics.txt \
	      summary/oxog_metrics.txt \
	      summary/hs_metrics.txt \
	      summary/gc_metrics.txt \
	      summary/gc_summary.txt
	      
TARGETS_LIST ?= $(HOME)/share/lib/resource_files/MSK-IMPACT-v4.sorted.list
	      
define idx-metrics
metrics/$1.idx_stats.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(BAM_INDEX) \
					    INPUT=$$(<) \
					    > $$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call idx-metrics,$(sample))))
					    
define aln-metrics
metrics/$1.aln_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(COLLECT_ALIGNMENT_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call aln-metrics,$(sample))))

define insert-metrics
metrics/$1.insert_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(COLLECT_INSERT_METRICS) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    HISTOGRAM_FILE=metrics/$1.insert_metrics.pdf \
					    MINIMUM_PCT=0.5")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call insert-metrics,$(sample))))

define oxog-metrics
metrics/$1.oxog_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(COLLECT_OXOG_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call oxog-metrics,$(sample))))

define hs-metrics
metrics/$1.hs_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(CALC_HS_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    BAIT_INTERVALS=$$(TARGETS_LIST) \
					    TARGET_INTERVALS=$$(TARGETS_LIST)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call hs-metrics,$(sample))))

define gc-metrics
metrics/$1.gc_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(COLLECT_GC_BIAS) \
					    INPUT=$$(<) \
					    OUTPUT=metrics/$1.gc_bias.txt \
					    CHART_OUTPUT=metrics/$1.gc_metrics.pdf \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    SUMMARY_OUTPUT=$$(@)")
					   
endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call gc-metrics,$(sample))))
		
summary/idx_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 1 --sample_names '$(SAMPLES)'")

summary/aln_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).aln_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 2 --sample_names '$(SAMPLES)'")

summary/insert_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 3 --sample_names '$(SAMPLES)'")

summary/oxog_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 4 --sample_names '$(SAMPLES)'")

summary/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 5 --sample_names '$(SAMPLES)'")

summary/gc_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 6 --sample_names '$(SAMPLES)'")
					  
summary/gc_summary.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_metrics.txt)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bam_metrics.R --option 7 --sample_names '$(SAMPLES)'")


..DUMMY := $(shell mkdir -p version; \
	     echo "picard" >> version/bam_metrics.txt; \
	     $(PICARD) CollectAlignmentSummaryMetrics --version &>> version/bam_metrics.txt; \
	     $(PICARD) CollectInsertSizeMetrics --version &>> version/bam_metrics.txt; \
	     $(PICARD) CollectOxoGMetrics --version &>> version/bam_metrics.txt; \
	     $(PICARD) CollectHsMetrics --version &>> version/bam_metrics.txt; \
	     $(PICARD) CollectGcBiasMetrics --version &>> version/bam_metrics.txt; \
             R --version >> version/bam_metrics.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: bam_metrics
