include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/deseq.$(NOW)

DESEQ_RNW = modules/rnaseq/deseq.Rnw
SWEAVE = $(RSCRIPT) modules/scripts/Sweave.R

CONDITION ?= condition
REF_CONDITION ?= ref

# pheno file: sample\tpheno with header
PHENO_FILE ?= pheno.txt

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all

deseq_results.txt : sumreads/geneCounts.txt
	mkdir -p graphics; $(SWEAVE) $(DESEQ_RNW) --condition $(CONDITION) --refCondition $(REF_CONDITION) --outFile $@ $< $(PHENO_FILE)


