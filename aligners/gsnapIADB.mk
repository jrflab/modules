include modules/Makefile.inc

OPTS = -d IADB -D ${GSNAP_REF} -B 4 -t 4 -A sam --novelsplicing=1 --pairexpect=200 -n 1 --quiet-if-excessive --nofails
ifeq ($(BAM_PHRED64),true)
	OPTS += -J 64 -j -31
endif
GSNAP_SGE_RREQ = $(call MEM_FREE,2G,4G) -q all.q -pe $(PARALLEL_ENV) 4 -now n

REQUIRED_FLAGS = 4
BAM_FILTER_FLAGS = 1536

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

VPATH = bam

LOGDIR = iadb/log

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach sample,$(SAMPLES),iadb/bam/$(sample).bam)

iadb/unaln_bam/%.bam : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" $(MKDIR) $(@D) $(LOGDIR); $(SAMTOOLS) view -f $(REQUIRED_FLAGS) -F $(BAM_FILTER_FLAGS) -bh $< > $@ 2> $(LOGDIR)/$(@F).log


iadb/fastq/%.1.fastq iadb/fastq/%.2.fastq : iadb/unaln_bam/%.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,5G,7G)" $(MKDIR) $(@D) $(LOGDIR); \
	$(call SAM_TO_FASTQ_MEM,5G) I=$< FASTQ=iadb/fastq/$*.1.fastq SECOND_END_FASTQ=iadb/fastq/$*.2.fastq &> $(LOGDIR)/$(@F).log

iadb/bam/%.gsnap.bam : iadb/fastq/%.1.fastq iadb/fastq/%.2.fastq
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1.5G,2G) -pe $(PARALLEL_ENV) 4" $(MKDIR) $(@D) $(LOGDIR); \
	$(MKDIR) $(@D) $(LOGDIR); $(GSNAP) $(OPTS) --read-group-id=$* $^ 2> $(LOGDIR)/$(@F).log | $(SAMTOOLS) view -bhS - > $@

iadb/bam/%.bam : iadb/bam/%.gsnap.sorted.filtered.markdup.bam
	$(MKDIR) $(@D); ln -v $< $@

include modules/bam_tools/processBam.mk
