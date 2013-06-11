# merge sample fastqs
include ~/share/modules/Makefile.inc

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

all : $(foreach i,1 2,$(foreach sample,$(SAMPLES),fastq/$(sample).$i.fastq.gz))

fastq/%.1.fastq.gz : fastq/%_*_R1_*.fastq.gz
	zcat $^ | gzip -c > $@ && $(RM) $^

fastq/%.2.fastq.gz : fastq/%_*_R2_*.fastq.gz
	zcat $^ | gzip -c > $@ && $(RM) $^
