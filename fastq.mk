# This module is used to extract fastq files from bam files

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

VPATH ?= unprocessed_bam

LOGDIR ?= log/fastq.$(NOW)

EXTRACT_TOOL ?= picard

ifeq ($(TRIM_READS),true)
   FASTQ_FILTER := trim
   TRIM_LENGTH ?= 150
   TRIM_OPTS ?= -l $(TRIM_LENGTH)
endif

FASTQ_TRIMMER = $(PERL) $(HOME)/share/scripts/trimFastq.pl

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq

ifdef UNSPLIT_SAMPLES
fastq: $(foreach split,$(UNSPLIT_SAMPLES),fastq/$(split).1.fastq.gz.md5)
else
fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz.md5)
endif

ifdef FASTQ_FILTER
fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5 : unprocessed_fastq/%.1.$(FASTQ_FILTER).fastq.gz.md5 unprocessed_fastq/%.2.$(FASTQ_FILTER).fastq.gz.md5
	$(INIT) ln $(<M) fastq/$*.1.fastq.gz && ln $(word 2,$(^M)) fastq/$*.2.fastq.gz && cp  $< fastq/$*.1.fastq.gz.md5 && cp $(word 2,$^) fastq/$*.2.fastq.gz.md5
else
fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5 : unprocessed_fastq/%.1.fastq.gz.md5 unprocessed_fastq/%.2.fastq.gz.md5
	$(INIT) ln $(<M) fastq/$*.1.fastq.gz && ln $(word 2,$(^M)) fastq/$*.2.fastq.gz && cp  $< fastq/$*.1.fastq.gz.md5 && cp $(word 2,$^) fastq/$*.2.fastq.gz.md5
endif

%.fastq.gz.md5 : %.fastq.gz
	$(call LSCRIPT,"$(MD5)")

ifeq (${EXTRACT_TOOL},picard)
fastq/%.1.fastq.gz.md5 fastq/%.2.fastq.gz.md5 : unprocessed_bam/%.bam
	$(call LSCRIPT_MEM,8G,10G,"TEMP=`mktemp`; mkfifo $${TEMP}_1; mkfifo $${TEMP}_2; \
	gzip < $${TEMP}_1 > fastq/$*.1.fastq.gz & \
	gzip < $${TEMP}_2 > fastq/$*.2.fastq.gz & \
	QUIET=true I=$< FASTQ=$${TEMP}_1 SECOND_END_FASTQ=$${TEMP}_2 && \
	md5sum fastq/$*.1.fastq.gz > fastq/$*.1.fastq.gz.md5 && \
	md5sum fastq/$*.2.fastq.gz > fastq/$*.2.fastq.gz.md5")
else
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : %.bam
	$(call LSCRIPT_MEM,10G,40G,"$(BAM2FASTQ) -o fastq/$*#.fastq $< &> $(LOG).bam2fastq.log && gzip < fastq/$*_1.fastq > fastq/$*.1.fastq.gz && gzip < fastq/$*_2.fastq > fastq/$*.2.fastq.gz")
endif

unprocessed_fastq/%.trim.fastq.gz.md5 : unprocessed_fastq/%.fastq.gz.md5
	$(call LSCRIPT_MEM,2G,3G,"$(CHECK_MD5) zcat $(<M) | $(FASTQ_TRIMMER) $(TRIM_OPTS) | gzip -c > $(@M) && $(MD5)")

unprocessed_fastq/%.readtrim.1.fastq.gz unprocessed_fastq/%.readtrim.2.fastq.gz : %.bam %.read_len
	$(call LSCRIPT_MEM,10G,12G,"NUM_READS=`awk '{ sum += $$1 } END { print sum }' $(word 2,$^)`; \
	MAX_LENGTH=`sort -k 2 $(word 2,$^) | awk -v nreads="$$NUM_READS" '$$1 / nreads > 0.4 { print $$2 }' | head -1`; \
	if [ "$$MAX_LENGTH" = "" ]; then MAX_LENGTH=`cut -d' ' -f 2 $(word 2,$^) | head -1`; fi; \
	TEMP=`mktemp`; mkfifo $${TEMP}_1; mkfifo $${TEMP}_2; \
	gzip < $${TEMP}_1 > fastq/$*.readtrim.1.fastq.gz & \
	gzip < $${TEMP}_2 > fastq/$*.readtrim.2.fastq.gz & \
	$(SAM_TO_FASTQ) I=$< FASTQ=$${TEMP}_1 SECOND_END_FASTQ=$${TEMP}_2 READ1_MAX_BASES_TO_WRITE=$$MAX_LENGTH READ2_MAX_BASES_TO_WRITE=$$MAX_LENGTH &> $(LOG)")

%.read_len : %.bam
	$(call LSCRIPT_MEM,4G,6G,"$(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' | awk ' > $@")

define merged-fastq
unprocessed_fastq/$1.%.fastq.gz.md5 : $$(foreach split,$2,unprocessed_fastq/$$(split).%.fastq.gz.md5)
	$$(INIT) $$(CHECK_MD5) zcat $^ | gzip > $@ 2> $$(LOG) && $$(MD5)
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-fastq,$(sample),$(split_lookup.$(sample)))))

