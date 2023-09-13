include modules/Makefile.inc

LOGDIR ?= log/wgs_metrics.$(NOW)

wgs_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).idx_stats.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).aln_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).insert_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).gc_metrics_summary.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics.txt) \
	      $(foreach sample,$(SAMPLES),metrics/$(sample).duplicate_metrics.txt) \
	      summary/idx_metrics.txt \
	      summary/aln_metrics.txt \
	      summary/insert_metrics.txt \
	      summary/oxog_metrics.txt \
	      summary/gc_metrics.txt \
	      summary/wgs_metrics.txt \
	      summary/duplicate_metrics.txt
	    
SAMTOOLS_THREADS = 4
SAMTOOLS_MEM_THREAD = 1G

GATK_THREADS = 4
GATK_MEM_THREAD = 2G

define picard-metrics
metrics/$1.idx_stats.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(BAM_INDEX) \
							INPUT=$$(<) \
							> $$(@)")
									   
metrics/$1.aln_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_ALIGNMENT_METRICS) \
							REFERENCE_SEQUENCE=$$(REF_FASTA) \
							INPUT=$$(<) \
							OUTPUT=$$(@)")
									   
metrics/$1.insert_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_INSERT_METRICS) \
							INPUT=$$(<) \
							OUTPUT=$$(@) \
							HISTOGRAM_FILE=metrics/$1.insert_metrics.pdf \
							MINIMUM_PCT=0.05")
									   
metrics/$1.oxog_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_OXOG_METRICS) \
							REFERENCE_SEQUENCE=$$(REF_FASTA) \
							INPUT=$$(<) \
							OUTPUT=$$(@)")
					    
metrics/$1.gc_metrics_summary.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_GC_BIAS) \
							INPUT=$$(<) \
							OUTPUT=metrics/$1.gc_metrics.txt \
							CHART_OUTPUT=metrics/$1.gc_metrics.pdf \
							REFERENCE_SEQUENCE=$$(REF_FASTA) \
							SUMMARY_OUTPUT=$$(@)")
					   
metrics/$1.wgs_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_WGS_METRICS) \
							INPUT=$$(<) \
							OUTPUT=$$(@) \
							REFERENCE_SEQUENCE=$$(REF_FASTA)")
							
metrics/$1.duplicate_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G -w 24:00:00,"set -o pipefail && \
							$$(COLLECT_DUP_METRICS) \
							INPUT=$$(<) \
							METRICS_FILE=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call picard-metrics,$(sample))))
	
summary/idx_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 1 --sample_names '$(SAMPLES)'")
					  
summary/aln_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).aln_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 2 --sample_names '$(SAMPLES)'")

summary/insert_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 3 --sample_names '$(SAMPLES)'")
					  
summary/oxog_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 4 --sample_names '$(SAMPLES)'")
					  
summary/gc_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_metrics_summary.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 5 --sample_names '$(SAMPLES)'")
					  
summary/wgs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 6 --sample_names '$(SAMPLES)'")
					  
summary/duplicate_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).duplicate_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/wgs_metrics.R --option 7 --sample_names '$(SAMPLES)'")

#..DUMMY := $(shell mkdir -p version; \
#	     $(SAMTOOLS) --version >> version/wgs_metrics.txt; \
#	     echo "gatk3" >> version/wgs_metrics.txt; \
#	     $(GATK) --version >> version/wgs_metrics.txt; \
#	     echo "picard" >> version/wgs_metrics.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: wgs_metrics
