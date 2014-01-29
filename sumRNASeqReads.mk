# This module deals with summarizing RNA-Seq reads over ensembl models and generating RPKM.  It requires that pre-processing of the bam files are done by GATK as input.
# Authors: Fong Chun Chan <fongchunchan@gmail.com> & Raymond Lim <raylim@mm.st>
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/sumReads.$(NOW)

SUM_READS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReads.R
SUM_EXONS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByExon.R
SUM_INTRONS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByIntron.R

HAS_CHR ?= false
ifeq ($(HAS_CHR),true)
SUM_READS_OPTS =
else
SUM_READS_OPTS = -r
endif


#SUM_INTRONS_BY_TRANSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByIntronByTranscript.R

#INTRONS_NON_OVERLAPPING_WITH_EXONS_FILE = ~/share/references/introns.non.overlapping.with.exons.txt

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all sumreads sumreads.intron.window sumexons sumintrons

all : sumreads # sumexons sumintrons

sumreads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.txt) sumreads/rpkm.txt sumreads/countsByGene.txt sumreads/countsByExon.txt
#sumexons : $(foreach sample,$(SAMPLES),sumexons/$(sample).sumexons.txt)
#sumintrons : $(foreach sample,$(SAMPLES),sumintrons/$(sample).sumintrons.txt)

sumreads/%.sumreads.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,15G,30G,"$(SUM_READS_RSCRIPT) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/rpkm.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/countsByGene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByGene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/countsByExon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByExon/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done



include ~/share/modules/processBamMD5.mk
