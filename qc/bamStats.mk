include modules/Makefile.inc
LOGDIR ?= log/bam_stats.$(NOW)

BAM_STATS_USE_REF ?= true

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: bam_stats

bam_stats: $(foreach sample,$(SAMPLES),metrics/$(sample).bam_stats.html)

metrics/%.bc : bam/%.bam
	$(call RUN,-s 6G -m 8G,"samtools stats $(if $(findstring true,$(BAM_STATS_USE_REF)),-r $(REF_FASTA)) $< > $@")
metrics/%.bam_stats.html : metrics/%.bc
	$(call RUN,-s 6G -m 8G,"plot-bamstats -p $(@D)/$* $<")
