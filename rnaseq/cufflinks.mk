# This module is used for running cufflinks
# input: $(SAMPLES) 
# Options: PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include ~/share/modules/Makefile.inc

LOGDIR = log/cufflinks.$(NOW)


NUM_CORES ?= 4
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks
CUFFCOMPARE = $(HOME)/share/usr/bin/cuffcompare
CUFFLINKS_OPTS = -b $(REF_FASTA) -g $(GENES_GTF) -p $(NUM_CORES) -u --no-update-check
CUFFCOMPARE_OPTS = -s $(REF_FASTA) -r $(GENES_GTF) -V


..DUMMY := $(shell mkdir -p version; $(CUFFLINKS) > version/tophat.txt; echo "options: $(CUFFLINKS_OPTS)" >> version/cufflinks.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all_cufflinks cufflinks cuffcmp

all_cufflinks : cufflinks cuffcmp
cufflinks : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
cuffcmp : cufflinks/cuffcmp/cc.stats

cufflinks/gtf/%.transcripts.gtf cufflinks/fpkm_tracking/%.isoforms.fpkm_tracking cufflinks/fpkm_tracking/%.genes.fpkm_tracking : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,$(NUM_CORES),2G,3G,"${CUFFLINKS} ${CUFFLINKS_OPTS} -o cufflinks/$* $< \
&& ln cufflinks/$*/transcripts.gtf cufflinks/gtf/$*.transcripts.gtf \
&& ln cufflinks/$*/isoforms.fpkm_tracking cufflinks/fpkm_tracking/$*.isoforms.fpkm_tracking \
&& ln cufflinks/$*/genes.fpkm_tracking cufflinks/fpkm_tracking/$*.genes.fpkm_tracking")

cufflinks/cuffcmp/cc.stats : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
	$(call LSCRIPT_MEM,10G,20G,"$(CUFFCOMPARE) $(CUFFCOMPARE_OPTS) -o $(@:.stats=) $^")
