# vim: set ft=make :
# pyloh loss of heterozygosity

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR := log/bwa.$(NOW)

PYLOH = PYTHONPATH=$(HOME)/usr/lib/python:$(HOME)/usr/lib/python2.7 $(HOME)/usr/bin/python $(HOME)/usr/bin/PyLOH.py

#BICSEQ = 



define pyloh-preprocess
pyloh/preprocess/$1_$2.PyLOH.counts : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,1G,2G,6,$$(PYLOH) preprocess $$(REF_FASTA) $$< $$(<<) $$(@D)/$1_$2 
endef
