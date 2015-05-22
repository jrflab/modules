# various bam processing steps
##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

LOGDIR ?= log/fixRG.$(NOW)

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
fixed_bams : $(addsuffix .md5,$(BAMS)) $(addsuffix .bai,$(BAMS))

bam/%.bam.md5 : unprocessed_bam/%.rg.bam.md5
	$(INIT) cp $< $@ && ln -f $(<M) $(@M)


include modules/bam_tools/processBam.mk
