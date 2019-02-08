include modules/Makefile.inc

LOGDIR ?= log/bam_to_fasta.$(NOW)
PHONY += unmapped_reads

bam_to_fasta : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).fasta)

define bam-to-fasta
unmapped_reads/%.fasta : unmapped_reads/%.bam
	$(call RUN,-n 4 -s 4G -m 9G,"$(SAMTOOLS2) fasta $$< > unmapped_reads/$$*.fasta")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fasta,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
