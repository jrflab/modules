include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

LOGDIR ?= log/fixRG.$(NOW)

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
fixed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : unprocessed_bam/%.rg.bam
	$(INIT) ln -f $(<) $(@)


include modules/bam_tools/processBam.mk
