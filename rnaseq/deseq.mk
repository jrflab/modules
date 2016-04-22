include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/deseq.$(NOW)

DESEQ_RNW = modules/rnaseq/deseq.Rnw
SWEAVE = $(RSCRIPT) modules/scripts/Sweave.R

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref

# pheno file: sample\tpheno with header
PHENO_FILE ?= pheno.txt

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all

deseq_results.txt : sumreads/geneCounts.txt
	mkdir -p graphics; $(SWEAVE) $(DESEQ_RNW) --condition $(DESEQ_CONDITION) --refCondition $(DESEQ_REF_CONDITION) --outFile $@ $< $(PHENO_FILE)


