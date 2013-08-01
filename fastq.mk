# This module is used to extract fastq files from bam files

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

VPATH ?= unprocessed_bam

LOGDIR ?= log/fastq.$(NOW)

EXTRACT_TOOL ?= picard

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq

fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)


ifeq (${EXTRACT_TOOL},picard)
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : unprocessed_bam/%.bam
	echo "TEMP=`mktemp`; mkfifo $${TEMP}_1; mkfifo $${TEMP}_2; \
	gzip < $${TEMP}_1 > fastq/$*.1.fastq.gz & \
	gzip < $${TEMP}_2 > fastq/$*.2.fastq.gz & \
	QUIET=true I=$< FASTQ=$${TEMP}_1 SECOND_END_FASTQ=$${TEMP}_2 &> $(LOGDIR)/$(@F).log" | $(LSCRIPT_MEM 8G)
else
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : %.bam
	$(call LAUNCH_MEM,10G,40G) $(BAM2FASTQ) -o fastq/$*#.fastq $< &> $(LOGDIR)/$(@F).bam2fastq.log && gzip < fastq/$*_1.fastq > fastq/$*.1.fastq.gz && gzip < fastq/$*_2.fastq > fastq/$*.2.fastq.gz
endif

fastq/%.readtrim.1.fastq.gz fastq/%.readtrim.2.fastq.gz : %.bam %.read_len
	echo "NUM_READS=`awk '{ sum += $$1 } END { print sum }' $(word 2,$^)`; \
	MAX_LENGTH=`sort -k 2 $(word 2,$^) | awk -v nreads="$$NUM_READS" '$$1 / nreads > 0.4 { print $$2 }' | head -1`; \
	if [ "$$MAX_LENGTH" = "" ]; then MAX_LENGTH=`cut -d' ' -f 2 $(word 2,$^) | head -1`; fi; \
	TEMP=`mktemp`; mkfifo $${TEMP}_1; mkfifo $${TEMP}_2; \
	gzip < $${TEMP}_1 > fastq/$*.readtrim.1.fastq.gz & \
	gzip < $${TEMP}_2 > fastq/$*.readtrim.2.fastq.gz & \
	$(SAM_TO_FASTQ) I=$< FASTQ=$${TEMP}_1 SECOND_END_FASTQ=$${TEMP}_2 READ1_MAX_BASES_TO_WRITE=$$MAX_LENGTH READ2_MAX_BASES_TO_WRITE=$$MAX_LENGTH &> $(LOGDIR)/$(@F).log" | $(LSCRIPT_MEM 10G,12G)

%.read_len : %.bam
	$(call LAUNCH_MEM,4G,5G) "$(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' | awk ' > $@"

