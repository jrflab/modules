# This module deals with summarizing RNA-Seq reads over ensembl models and generating RPKM.  It requires that pre-processing of the bam files are done by GATK as input.
# Authors: Fong Chun Chan <fongchunchan@gmail.com> & Raymond Lim <raylim@mm.st>
include ~/share/modules/Makefile.inc
include ~/share/modules/hg19.inc

SAMPLE_FILE ?= samples.txt
LOGDIR ?= log

SAMPLES = $(shell cat $(SAMPLE_FILE))
SUM_READS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReads.R
SUM_EXONS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByExon.R
SUM_INTRONS_RSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByIntron.R
#SUM_INTRONS_BY_TRANSCRIPT = ${RSCRIPT} ~/share/scripts/summarizeRNASeqReadsByIntronByTranscript.R

INTRONS_NON_OVERLAPPING_WITH_EXONS_FILE = ~/share/references/introns.non.overlapping.with.exons.txt

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all sumreads sumreads.intron.window sumexons sumintrons

INTRON_WINDOW = 50

all : sumreads sumreads.intron.window sumexons sumintrons

sumreads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.txt)
sumreads_introns_non_overlapping_with_exons : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.introns.non.overlapping.with.exons.txt)
sumreads.intron.window : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.intron.window.$(INTRON_WINDOW).txt)

sumexons : $(foreach sample,$(SAMPLES),sumexons/$(sample).sumexons.txt)

sumintrons : $(foreach sample,$(SAMPLES),sumintrons/$(sample).sumintrons.txt)
sumintrons_introns_non_overlapping_with_exons : $(foreach sample,$(SAMPLES),sumintrons/$(sample)_sumintrons_non_overlapping_with_exons.txt)
sumintrons.intron.window : $(foreach sample,$(SAMPLES),sumintrons/$(sample).sumintrons.intron.window.$(INTRON_WINDOW).txt)

####
# Get the gene-centric summarized reads
####
define getSumReads 
$1 : $2 $3
	$$(call INIT_MEM,15G,30G) $$(SUM_READS_RSCRIPT) $4 $$(TXDB) $$< $$@ &> $$(LOG)
endef

$(foreach sample,$(SAMPLES),$(eval $(call getSumReads,sumreads/$(sample).sumreads.txt,bam/$(sample).bam,bam/$(sample).bam.bai,$(SUM_READS_FLAGS))))
$(foreach sample,$(SAMPLES),$(eval $(call getSumReads,sumreads/$(sample).sumreads.intron.window.$(INTRON_WINDOW).txt,bam/$(sample).bam,bam/$(sample).bam.bai,-a -w $(INTRON_WINDOW))))
$(foreach sample,$(SAMPLES),$(eval $(call getSumReads,sumreads/$(sample).sumreads.introns.non.overlapping.with.exons.txt,bam/$(sample).bam,bam/$(sample).bam.bai,-a -i $(INTRONS_NON_OVERLAPPING_WITH_EXONS_FILE))))

####
# Get the intronic-centric summarized reads
####
define getSumIntrons
$1 : $2 $3
	SGE_RREQ="$$(SGE_RREQ) $$(call MEM_FREE,15G,30G)" \
	mkdir -p $$(@D)/logs;\
	$$(SUM_INTRONS_RSCRIPT) $4 $$(TXDB) $$< $$@ &> $$(@D)/logs/$$(@F).log
endef

$(foreach sample,$(SAMPLES),$(eval $(call getSumIntrons,sumintrons/$(sample).sumintrons.txt,bam/$(sample).bam,bam/$(sample).bam.bai,-a)))
$(foreach sample,$(SAMPLES),$(eval $(call getSumIntrons,sumintrons/$(sample)_sumintrons_non_overlapping_with_exons.txt,bam/$(sample).bam,bam/$(sample).bam.bai,-a -i $(INTRONS_NON_OVERLAPPING_WITH_EXONS_FILE))))

####
# Get the exon-centric summarized reads
####
define getSumExons
$1 : $2 $3
	$$(call INIT_MEM,15G,30G) $$(SUM_EXONS_RSCRIPT) $4 $$(TXDB) $$< $$@ &> $$(LOG)
endef

$(foreach sample,$(SAMPLES),$(eval $(call getSumExons,sumexons/$(sample).sumexons.txt,bam/$(sample).bam,bam/$(sample).bam.bai,$(SUM_EXONS_FLAGS))))

include ~/share/modules/processBam.mk
