# merge sample fastqs
include ~/share/modules/Makefile.inc

LOGDIR = log/merge_fastq.$(NOW)

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:
.SECONDEXPANSION: 

all : $(foreach i,1 2,$(foreach sample,$(SAMPLES),fastq/$(sample).$i.fastq.gz.md5))

fastq/%.1.fastq.gz.md5 : fastq/$$*_*_R1_*.fastq.gz
	$(call LSCRIPT,"zcat $(sort $^) | gzip -c > $(@:.md5=) && $(MD5) && $(RM) $^")

#$(INIT) echo "$^" | tr " " "\n" | sort | xargs zcat | gzip -c > $@ && $(RM) $^

fastq/%.2.fastq.gz.md5 : fastq/$$*_*_R2_*.fastq.gz
	$(call LSCRIPT,"zcat $(sort $^) | gzip -c > $(@:.md5=) && $(MD5) && $(RM) $^")
