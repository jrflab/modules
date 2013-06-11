# vim: set ft=make :
# Perform rename operations on miseq fastq files

include ~/share/modules/Makefile.inc

SAMPLES = $(notdir $(subst _,,$(patsubst %_L001_R1_001.fastq.gz,%,$(wildcard fastq/*_L001_R1_001.fastq.gz))))

.PHONY : fastq all

all : fastq samples.txt

fastq :
	cd fastq; rename _L001_R1_001.fastq.gz .1.fastq.gz *.fastq.gz; rename _L001_R2_001.fastq.gz .2.fastq.gz *.fastq.gz; rename _S S *.fastq.gz; cd ..; \
	SMALL=`find fastq/* -size -1000 | sed 's/\..*.fastq.gz//'`; \
	for x in $$SMALL; do $(RM) $x.*.fastq.gz; done 


samples.txt : fastq
	ls fastq/*.1.fastq.gz | sed 's;.*/;;; s/\..*//' > $@
