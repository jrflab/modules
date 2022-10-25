include modules/Makefile.inc

LOGDIR = log/sum_reads.$(NOW)

sumreads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_gene.txt) \
	   $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_exon.txt) \
	   sumreads/rpkm_by_gene.txt \
	   sumreads/rpkm_by_exon.txt \
	   sumreads/counts_by_gene.txt \
	   sumreads/counts_by_exon.txt

SUM_READS_OPTS =
REF ?= b37

sumreads/%.sumreads.by_gene.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-v $(SUMREADS_ENV) -s 24G -m 48G,"$(SUM_READS_RSCRIPT) --genome $(REF) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/%.sumreads.by_exon.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-v $(SUMREADS_ENV) -s 24G -m 48G,"$(SUM_EXONS_RSCRIPT) --genome $(REF) --outFile $@ $(SUM_READS_OPTS) $<")

sumreads/rpkm_by_gene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_gene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/rpkm_by_exon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_exon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 6 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_by_gene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_gene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByGene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_by_exon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads.by_exon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/exonCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

..DUMMY := $(shell mkdir -p version; \
	     $(SUMREADS_ENV)/bin/R --version >> version/sumreads.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: sumreads
