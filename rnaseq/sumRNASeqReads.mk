include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/sum_reads.$(NOW)

DEFAULT_ENV = $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.6

SUM_READS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReads.R
SUM_EXONS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReadsByExon.R
SUM_INTRONS_RSCRIPT = ${RSCRIPT} modules/rnaseq/summarizeRNASeqReadsByIntron.R
SUM_READS_OPTS =

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all sumreads

SUM_TYPE = byGene byExon

all : $(foreach type,$(SUM_TYPE),$(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.$(type).txt)) sumreads/rpkm_by_gene.txt sumreads/rpkm_by_exon.txt sumreads/counts_by_gene.txt sumreads/counts_by_exon.txt

sumreads/%.sumreads.byGene.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-v $(DEFAULT_ENV) -s 24G -m 48G,"$(SUM_READS_RSCRIPT) --genome $(REF) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/%.sumreads.byExon.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-v $(DEFAULT_ENV) -s 24G -m 48G,"$(SUM_EXONS_RSCRIPT) --genome $(REF) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/rpkm_by_gene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byGene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/rpkm_by_exon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byExon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 6 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_by_gene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byGene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByGene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_by_exon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.byExon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/exonCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

include modules/bam_tools/processBam.mk
