# This module deals with summarizing RNA-Seq reads over ensembl models and generating RPKM.  It requires that pre-processing of the bam files are done by GATK as input.
# Authors: Fong Chun Chan <fongchunchan@gmail.com> & Raymond Lim <raylim@mm.st>
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/sumReads.$(NOW)

SUM_READS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReads.R
SUM_EXONS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReadsByExon.R
SUM_INTRONS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReadsByIntron.R

SUM_READS_OPTS =


#SUM_INTRONS_BY_TRANSCRIPT = ${RSCRIPT} scripts/summarizeRNASeqReadsByIntronByTranscript.R

#INTRONS_NON_OVERLAPPING_WITH_EXONS_FILE = ~/share/references/introns.non.overlapping.with.exons.txt

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all sumreads

SUM_TYPE = byGene byExon

all : $(foreach type,$(SUM_TYPE),$(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.$(type).txt)) sumreads/geneRPKM.txt sumreads/exonRPKM.txt sumreads/geneCounts.txt sumreads/exonCounts.txt

#sumintrons : $(foreach sample,$(SAMPLES),sumintrons/$(sample).sumintrons.txt)

sumreads/%.sumreads.byGene.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,30G,60G,"$(SUM_READS_RSCRIPT) -g $(REF) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/%.sumreads.byExon.txt : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,20G,40G,"$(SUM_EXONS_RSCRIPT) --txdb $(ENSEMBL_TXDB) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/geneRPKM.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byGene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/exonRPKM.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byExon.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 6 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/geneCounts.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byGene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByGene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/exonCounts.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byExon.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/exonCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done


include modules/bam_tools/processBam.mk
