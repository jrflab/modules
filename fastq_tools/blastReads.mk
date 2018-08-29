include modules/Makefile.inc

LOGDIR ?= log/blast_reads.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: blast_reads

VPATH = unmapped_reads

blast_reads : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).blast)

unmapped_reads/%.blast : unmapped_reads/%.fasta
	$(call RUN,-n 32 -s 3G -m 4G -w 360,"blastn -num_threads 32 -evalue 0.001 -word_size 28 -db ~/share/reference/ncbi_nt/nt -query $< -outfmt 7 -out unmapped_reads/$*.blast")
