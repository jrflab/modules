include modules/Makefile.inc

LOGDIR ?= log/extract_unmapped.$(NOW)
PHONY += unmapped_reads

extract_unmapped : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).bam)

unmapped_reads/%.bam : %.bam
	$(call RUN,-n 4 -s 4G -m 9G,"$(SAMTOOLS2) view -f 0x04 -h -@ 4 -b $< -o unmapped_reads/$*.bam")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)