# This module is used for running cufflinks
# input: $(SAMPLES) 
# Options: PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include ~/share/modules/Makefile.inc

NUM_CORES ?= 4
CUFFLINKS_OPTS = -b $(REF_FASTA) -g $(REFSEQ_GTF) -p $(NUM_CORES) -u --no-update-check

SAMPLE_FILE ?= samples.txt

REF_GTF = $(HOME)/share/references/refseq.gtf

NUM_CORES = 4
CUFFLINKS_OPTS = -b $(REF_FA) -g $(REF_GTF) -p $(NUM_CORES) -u --no-update-check

SAMPLES = $(shell cat $(SAMPLE_FILE))

.PHONY : cufflinks cuffcmp

cufflinks : $(foreach sample,$(SAMPLES),cufflinks/$(sample).cufflinks.timestamp)
cuffcmp : cufflinks/cuffcmp.timestamp

cufflinks/%.cufflinks.timestamp : bam/%.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,5G,6G) -pe $(PARALLEL_ENV) $(NUM_CORES) -q all.q" \
	$(MKDIR) $(@D)/logs;\
	${CUFFLINKS} ${CUFFLINKS_OPTS} -o $(@D)/$* $< &> $(@D)/logs/$*.log && touch $@

cufflinks/cuffcmp.timestamp : $(foreach sample,$(SAMPLES),cufflinks/$(sample).cufflinks.timestamp)
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,20G) -q all.q" \
	$(MKDIR) $(@D)/logs;\
	cd $(@D);\
	echo -e "transFragID\tlocusID\trefGeneID\tclassCode" $(patsubst %,\\t%,$(SAMPLES)) > cuffcmp_header.txt;\
	$(CUFFCMP) -s $(REF_FASTA) -r $(REFSEQ_GTF) -V $(patsubst $(@D)/%.cufflinks.timestamp,%/transcripts.gtf,$^) &> logs/cuffcmp.log && touch $(@F)
