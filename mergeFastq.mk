# merge sample fastqs
include ~/share/modules/Makefile.inc

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))
LOGDIR = log/merge_fastq.$(NOW)

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:
.SECONDEXPANSION: 

all : $(foreach i,1 2,$(foreach sample,$(SAMPLES),fastq/$(sample).$i.fastq.gz))

fastq/%.1.fastq.gz : fastq/$$*_*_R1_*.fastq.gz
	$(INIT) echo "$^" | tr " " "\n" | sort | xargs zcat | gzip -c > $@ && $(RM) $^

fastq/%.2.fastq.gz : fastq/$$*_*_R2_*.fastq.gz
	$(INIT) echo "$^" | tr " " "\n" | sort |  xargs zcat | gzip -c > $@ && $(RM) $^
