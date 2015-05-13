# integrate only on rnaseq
# pre-req: tophat

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/integrate_rnaseq.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(INTEGRATE) &> version/integrate.txt")

INTEGRATE = $(HOME)/share/usr/bin/Integrate

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: integrate_rnaseq

integrate_rnaseq : $(foreach sample,$(SAMPLES),integrate/breakpoints/$(sample).breakpoints.tsv)

integrate/reads/%.reads.txt integrate/sum/%.sum.tsv integrate/exons/%.exons.tsv integrate/breakpoints/%.breakpoints.tsv : bam/%.bam.md5 bam/%.bam.bai
	$(call LSCRIPT_MEM,8G,60G,"mkdir -p integrate/reads integrate/sum integrate/exons integrate/breakpoints; $(INTEGRATE) fusion -reads integrate/reads/$*.reads.txt -sum integrate/sum/$*.sum.tsv -ex integrate/exons/$*.exons.tsv -bk integrate/breakpoints/$*.breakpoints.tsv $(REF_FASTA) $(INTEGRATE_ANN) $(INTEGRATE_BWTS) $(<M) $(<M)")

include modules/bam_tools/processBam.mk
